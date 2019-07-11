################################################################################
# This script generates an ASKAP spectral line cube QC HTML report. 
# This version runs on the final continuum subtracted cube/mosaic. 
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- slurmOutput/*sh
# 
#
# To run type:
#
# "python <script name> <main directory name> <SBID>
# e.g. python val_v1.0.py Eridanus_processed 1234
#
# Author: Bi-Qing For
# Email: biqing.for@icrar.org
# 
# Modified Date: 2 July 2019
#
################################################################################

import os.path
import sys
import glob
from astropy.io import fits
import numpy as np
from datetime import datetime
from PIL import Image
from astropy.stats import median_absolute_deviation
import math
from scipy.constants import k as k_B
import subprocess 

PI = math.pi

main_dir = sys.argv[1]
sbid = str(sys.argv[2])
askapsoft = '0.11.0'  # This is a random number until finding out where this information is provided
bmaj = '30' # This is a random number 
bmin = '30' # This is a random number

### Need to do checking if files exist

#################################################
#       Functions for the main program
#################################################

#def get_FitsHeader(fitsimage):
#    """
#    Getting basic information of processed mosaic
#    """

#    hdu = fits.getheader(fitsimage)
#    bmaj = round(float(hdu.header['BMAJ'])*3600., 2)  #arcsec
#    bmin = round(float(hdu.header['BMIN'])*3600., 2)  #arcsec
    
def get_Metadata(metafile):
    """
    Getting basic information on observed field (one field only).
    """

    mslist_file = open(metafile, 'r')
    LINES = mslist_file.readlines()
    mslist_file.close()

    obs_date = 'Observed from'
    code = 'Code'
    duration = 'Total elapsed time'
    antenna = 'antennas'

    for i in range(len(LINES)):
        line = LINES[i]
        if line.find(antenna) >=0:
            TOKS = line.split()
            n_ant = TOKS[5][-2:]
        if line.find(obs_date) >=0:
            TOKS = line.split()
            start_obs_date = TOKS[6]
            end_obs_date = TOKS[8]
        if line.find(duration) >=0:
            TOKS = line.split()
            tobs = float(TOKS[10]) # in second
        if line.find(code) >= 0:
            next_line = LINES[i+1]
            TOKS = next_line.split()
            field = TOKS[5]
            ra = TOKS[6][:-5]
            dec = TOKS[7][:-4]
            
    return n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec

def get_Metadata_freq(metafile_science):
    """
    Getting basic information on observed field (one field only). Frequency and channel related items. 
    """

    mslist_file = open(metafile_science, 'r')
    LINES = mslist_file.readlines()
    mslist_file.close()

    frame = 'Frame'
    
    for i in range(len(LINES)):
        line = LINES[i]
        if line.find(frame) >=0:
            next_line = LINES[i+1]
            TOKS = next_line.split()
            chan_width = float(TOKS[10])*1000. # convert kHz to Hz
            bw = TOKS[11]  #kHz
            cfreq = TOKS[12] #MHz

    return chan_width, bw, cfreq    
    

def get_Frequency_Range(cubestat_contsub):
    """
    Frequency range of the mosaic.
    """
    
    line = subprocess.check_output(['sed', '-n', '3p', cubestat_contsub])
    TOKS = line.split()
    start_freq = round(float(TOKS[1]), 3)
    
    line = subprocess.check_output(['tail', '-1', cubestat_contsub])
    TOKS = line.split()
    end_freq = round(float(TOKS[1]), 3)

    return start_freq, end_freq

def make_Thumbnail(image, thumb_img, sizeX, sizeY):
    """
    Making thumbnail image.
    """
    im = Image.open(image)
    im.thumbnail(sizeX, sizeY)
    im.save(main_dir+ '/' + thumb_img)

    return thumb_img

def cal_beam_MADMFD(infile):
    """
    Calculating the MAD of max flux density of each beam.
    """

    data = np.loadtxt(infile)
    maxfdensity = data[:,8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 2)

    return mad_maxfdensity
    
def qc_Max_Flux_Density (infile, delta_freq_range, mean_beamMADMFD):
    """
    Evaluating the max flux density (MFD) of the mosaic.
    Metric to check for the effect of solar interference. 
    First, check the processed mosaic contains at least 5 MHz of data. 
    Then, calculating the MAD of MFD and comparing it to the mean of MADMFD of central 16 beams. 
    """

    if delta_freq_range > 5.0:  # need to have at least 5 MHz bandwidth of data to check for meaningful variation
        data = np.loadtxt(infile)
        maxfdensity = data[:,8]
        mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 2)
        if (mad_maxfdensity > mean_beamMADMFD): # need to compare to the centre 16 beams
            maxfden_id = 'bad' 
        else:
            maxfden_id = 'good' 

    else:
        print 'Warning: Processed data are less than 5 MHz, variation of max flux density is not meaningful.'
        maxfden_id = 'uncertain'
        
    return mad_maxfdensity, maxfden_id


def cal_Theoretical_RMS (n_ant, tobs, chan_width):
    """
    Calculating the theoretical rms noise for ASKAP. Assuming natural weighting and not taking into account fraction of flagged data. 
    """

    tsys = 50       # K
    antdiam = 12    # m
    aper_eff = 0.8  # aperture efficiency
    coreff = 0.8    # correlator efficiency
    npol = 2.0      # Number of polarisation, npol = 2 for images in Stokes I, Q, U, or V
    
    anteff = PI*(antdiam/2)**2. * aper_eff
    SEFD = 2. * k_B * 1e26 *tsys/anteff   
    rms_jy = SEFD/(coreff*math.sqrt(npol*n_ant*(n_ant-1)*chan_width*tobs))

    return rms_jy


def get_Flagging (flagging_file):
    """
    Getting flagging statistics.
    """

    line = subprocess.check_output(['tail', '-1', flagging_file])
    TOKS = line.split()
    flagged_stat = TOKS[-1]

    return flagged_stat

"""
def check_Cleaning ():
    
    Checking if cleaning has performed down to the defined threshold. Both major and minor cycles. 
    
    This will involve ncore specified. 
"""    


###########################################################

# Required files 

metafile = sorted(glob.glob('metadata/mslist-*txt'))[0]
metafile_science = sorted(glob.glob('metadata/mslist-scienceData*txt'))[0]
cubestat_contsub = glob.glob(main_dir + '/cubeStats*contsub.txt')[0]
flagging_file = glob.glob('slurmOutput/flag.out.txt')[0]  # temporary file from Wasim


#############################    
# HTML related tasks
#############################

html_name = 'spectral_report_SB' + sbid + '.html'

### Making thumbnail images

sizeX = str(70)
sizeY = str(70)

cube_plots = glob.glob(main_dir + '/cubePlot*.png')
beamNoise_plots = glob.glob(main_dir + '/beamNoise*.png')
beamMinMax_plots = glob.glob(main_dir + '/beamMinMax*.png')

thumb_cubeplots = []
thumb_beamNoise = []
thumb_beamMinMax = []

for image in cube_plots:
    thumb_img = 'thumb_'+ sizeX + '_'+ image.partition('/')[2]
    thumb_cubeplots.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY)

for image in beamNoise_plots:
    thumb_img = 'thumb_'+ sizeX + '_'+ image.partition('/')[2]
    thumb_beamNoise.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY)

for image in beamMinMax_plots:
    thumb_img = 'thumb_'+ sizeX + '_'+ image.partition('/')[2]
    thumb_beamMinMax.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY)


# Calling the functions

n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec = get_Metadata(metafile)
chan_width, bw, cfreq = get_Metadata_freq(metafile_science)
start_freq, end_freq = get_Frequency_Range(cubestat_contsub)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.

tobs_hr = round(tobs/3600.,2) # convert tobs from second to hr
chan_width_kHz = round(chan_width/1000.,3) # convert Hz to kHz

freq_range = str(start_freq)+'--'+str(end_freq)
delta_freq_range = end_freq - start_freq

### Validation parameters

flagged_stat = get_Flagging(flagging_file)

beamstat_contsub = glob.glob(main_dir + '/' + field +'/cubeStats*beam??.contsub.txt')

for beamfile in beamstat_contsub:
    beam_num = beamfile[-14:-12]
    if int(beam_num) < 16:
        beamMADMFD = cal_beam_MADMFD(beamfile)

mean_beamMADMFD = np.mean(beamMADMFD)

QC_mad_maxfden, QC_maxfden_id = qc_Max_Flux_Density(cubestat_contsub, delta_freq_range, mean_beamMADMFD) #Continuum subtracted
#QC_maxfden_id = 'bad'
#QC_mad_maxfden = 2.4


##############################
# Creating html report
##############################

css_style = """<style>
                  body {
                        background-color: white;
                  }
                  h1 {
                      color: black;
                  }
                  p {
                     color: black;
                  }
                  th, td {
                         padding: 5px;
                  }
                  table {
                         width: "100%";
                         border: "1";
                         border-spacing: 2px;
                  }
                  img {
                       border: 1px solid #ddd; /* Gray border */
                       border-radius: 4px;
                       padding: 5px;
                  }
                  #good {
                         background-color:#00FF00;
                  }

                  #uncertain {
                         background-color:#FFD700;
                  }
                  #bad {
                         background-color:#CD5C5C;
                  }
                  </style>"""

html = open(html_name, 'w')
html.write("""<!DOCTYPE HTML>
              <html>
                    <head>
                    %s
                    </head>
                    <body>
                    <center>
                    <title>ASKAP WALLABY Spectral Line Data Validation Report</title>
                    <h1 align="middle">WALLABY Spectral Line Data Validation Report</h1>
                    """ % (css_style))

html.write("""<h2 align="middle">Observation</h2>
                    <table width="100%" border="1">
                    <tr>
                        <th>SBID</th>
                        <th>No. of Antennas</th>
                        <th>Obs Start Date/Time</th>
                        <th>Obs End Date/Time</th>   
                        <th>Duration<br>(hr)</th>
                        <th>Field</th>
                        <th>R.A.</th>
                        <th>Decl.</th>
                        <th>Total Bandwidth <br>(kHz)</th>
                        <th>Central Frequency<br>(MHz)</th>
                    </tr>
                    <tr align="middle">
                        <td>{0}</td>
                        <td>{1}</td>
                        <td>{2}</td>
                        <td>{3}</td>
                        <td>{4:.1f}</td>
                        <td>{5}</td>
                        <td>{6}</td>
                        <td>{7}</td>
                        <td>{8}</td>
                        <td>{9}""".format(sbid,
                                          n_ant,
                                          start_obs_date,
                                          end_obs_date,
                                          tobs_hr,
                                          field,
                                          ra,
                                          dec,
                                          bw,
                                          cfreq))

html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Processed Image Cube</h2>
                    <table width="100%" border="1">
                    <tr>
                        <th>ASKAPsoft<br>version</th>
                        <th>Frequency Range<br>(MHz)</th>
                        <th>Channel Width<br>(kHz)</th>
                        <th>Synthesised Beam bmaj<br>(arcsec)</th>
                        <th>Synthesised Beam bmin<br>(arcsec)</th>
                    </tr>
                    <tr align="middle">
                        <td>{0}</td>
                        <td>{1}</td>
                        <td>{2}</td>
                        <td>{3}</td>
                        <td>{4}""".format(askapsoft,
                                          freq_range,
                                          chan_width_kHz,
                                          bmaj,
                                          bmin))

html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Mosaic Statistics</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Image Cube</th>          
                    <th>Continuum Subtracted Cube</th>
                    <th>Residual Cube</th>
                    <th>Total Flagged Visibilities</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    </td>
                    <td>{12}
                    """.format(cube_plots[1],
                                   main_dir+'/'+ thumb_cubeplots[1],
                                   sizeX,
                                   sizeY,
                                   cube_plots[0],
                                   main_dir+'/'+ thumb_cubeplots[0],
                                   sizeX,
                                   sizeY,
                                   cube_plots[2],
                                   main_dir+'/'+ thumb_cubeplots[2],
                                   sizeX,
                                   sizeY,
                                   flagged_stat))


html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Beams Statistics</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Beam Image Cube</th>          
                    <th>Continuum Subtracted Beam Cube</th>
                    <th>Residual Beam Cube</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile</br></p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile</br></p>
                    </td>
                    <td>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile</br></p>
                    </td>
                    <tr align="middle">
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</br></p>
                    </td>
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</br></p>
                    </td>
                    <td>
                    <a href="{20}" target="_blank"><img src="{21}" width="{22}" height="{23}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</br></p>
                    """.format(beamMinMax_plots[1],
                                    main_dir+'/'+ thumb_beamMinMax[1],
                                    sizeX,
                                    sizeY,
                                    beamMinMax_plots[0],
                                    main_dir+'/'+ thumb_beamMinMax[0],
                                    sizeX,
                                    sizeY,
                                    beamMinMax_plots[2],
                                    main_dir+'/'+ thumb_beamMinMax[2],
                                    sizeX,
                                    sizeY,
                                    beamNoise_plots[0],
                                    main_dir+'/'+ thumb_beamNoise[0],
                                    sizeX,
                                    sizeY,
                                    beamNoise_plots[1],
                                    main_dir+'/'+ thumb_beamNoise[1],
                                    sizeX,
                                    sizeY,
                                    beamNoise_plots[2],
                                    main_dir+'/'+ thumb_beamNoise[2],
                                    sizeX,
                                    sizeY))

### Finish HTML report with generated time stamp

html.write("""</td> 
              </tr>
              </table>
              <h2 align="middle">Validation Metrics</h2>
              <table width="100%" border="1">
              <tr>
              <th>Max Flux Density <br>(Jy)</th>
              <th>item</th>
              <th>item</th>   
              </tr>
              <tr align="middle">
              <td id='{0}'>{1}</td>
              <td id='{2}'>{3}</td>
              <td id='{4}'>{5}""".format(QC_maxfden_id,
                                         QC_mad_maxfden,
                                         QC_maxfden_id,
                                         QC_mad_maxfden,
                                         QC_maxfden_id,
                                         QC_mad_maxfden))

html.write("""</td>
              </tr>
              </table>
                                             
              <p align=left>Generated at {0} <br>
              <i> Report bugs to 
              <a href="mailto:biqing.for@icrar.org">Bi-Qing For</a></i>
              </p>

              </body>
              </html>""".format(str(datetime.now())))





html.close()
print "Spectral line validation report written to '{0}'.".format(html_name)
                    
