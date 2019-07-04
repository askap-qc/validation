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
    frame = 'Frame'
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
        if line.find(frame) >=0:
            next_line = LINES[i+1]
            TOKS = next_line.split()
            chan0 = float(TOKS[9])
            chan_width = float(TOKS[10])*1000. # convert kHz to Hz
            bw = TOKS[11]
            cfreq = TOKS[12]
            
    return n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, chan0, chan_width, bw, cfreq 

## Might have to implement a statement to check if commented out # is used in the param file

def get_Channel_Range(param):
    """
    Checking if channel range is used for processing.
    """
    
    channel_range = 'CHAN_RANGE_SCIENCE'

    with open(param) as infile:
        for line in infile:
            if line.find(channel_range) >=0:
                chan_range = line.partition('=')[2]
                start_chan_range = chan_range.partition('-')[0]
                end_chan_range = chan_range.partition('-')[2]
                break
            else:
                chan_range = "Full bandwidth"
                start_chan_range = False
                end_chan_range = False

    return chan_range, start_chan_range, end_chan_range

def convert_Channel_Range(chan0, start_chan_range, end_chan_range, chan_width):
    """
    Converting channel range in channel to frequency.
    """
    
    freq_start_chan_range = round(chan0 + float(start_chan_range)*chan_width/1e6, 3) # MHz
    freq_end_chan_range = round(chan0 + float(end_chan_range)*chan_width/1e6, 3) # MHz

    return freq_start_chan_range, freq_end_chan_range


def make_Thumbnail(image, thumb_img, sizeX, sizeY):
    """
    Making thumbnail image.
    """
    im = Image.open(image)
    im.thumbnail(sizeX, sizeY)
    im.save(main_dir+ '/' + thumb_img)

    return thumb_img

def qc_Max_Flux_Density (infile, pro_freq_range):
    """
    Evaluate the max flux density.
    """
#    if pro_freq_range > 0:
        
    data = np.loadtxt(infile)
    rms, maxfdensity = data[:,3], data[:,8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 2)
#    med_rms = round(np.median(rms),2)
    
    if (mad_maxfdensity >= 2.0):
        maxfden_id = 'bad' 
    elif (mad_maxfdensity < 2.0):
        maxfden_id = 'good' 
    else:
        maxfden_id = 'uncertain' ### currently don't have a condition that's uncertain

#    if (med_rms <= expected_rms):
#        med_rms_id = 'good'
#    else:
#        med_rms_id = 'bad' ### might need to have one more condition for acceptable rms if > expected rms

    return mad_maxfdensity, maxfden_id


def cal_Theoretical_RMS (n_ant, tobs, chan_width):
    """
    Calculate theoretical rms noise for ASKAP. Assuming natural weighting and not taking into account fraction of flagged data. 
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

###########################################################

# Required files 

metafile = sorted(glob.glob('metadata/mslist-scienceData*txt'))[0]
param_file = sorted(glob.glob('slurmOutput/*.sh'))
cubestat_cont = glob.glob(main_dir+'/cubeStats*contsub.txt')[0]

# Check if there is on parameter input file exists in the slurmOutput directory.
# If it does, select the latest one. 

if param_file >=1:
    index = len(param_file)
    param = param_file[index-1]
else:
    param = param_file[0]

#############################    
# HTML related tasks
#############################

html_name = 'spectral_QC_report_SB' + sbid + '.html'

### Making thumbnail images

sizeX = str(120)
sizeY = str(120)

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

n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, chan0, chan_width, bw, cfreq = get_Metadata(metafile)
chan_range, start_chan_range, end_chan_range = get_Channel_Range(param)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.

print theoretical_rms_mjy

tobs_hr = round(tobs/3600.,2) # convert tobs from second to hr
chan_width_kHz = round(chan_width/1000.,3) # convert Hz to kHz

pro_freq_range = 0 # in the case of full bandwidth is used for processing
if chan_range != "Full bandwidth":
    freq_start_chan_range, freq_end_chan_range = convert_Channel_Range(chan0, start_chan_range, end_chan_range, chan_width)
    pro_freq_range = freq_end_chan_range - freq_start_chan_range
    chan_range = str(freq_start_chan_range)+'-'+str(freq_end_chan_range)

### Validation parameters

QC_mad_maxfden, QC_maxfden_id = qc_Max_Flux_Density(cubestat_cont, pro_freq_range) #Continuum subtracted



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
                        <td>{3:.1f}</td>
                        <td>{4}</td>
                        <td>{5}</td>
                        <td>{6}</td>
                        <td>{7}</td>
                        <td>{8}""".format(sbid,
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
                        <th>Channel Range<br>(MHz)</th>
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
                                          chan_range,
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
                                   sizeY))


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
                    
