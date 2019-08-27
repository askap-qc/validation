################################################################################
# This script generates an ASKAP spectral line cube HTML report. 
# This version runs on the final continuum subtracted cube/mosaic. 
#
# Compatibility: Python version 3.x
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- slurmOutput/*sh
# -- main_dir/image.restored.i.SB<SBID>.cube.contsub.fits
#
# To run type:
#
# "python <script name> <main directory name> <SBID>
# e.g. python ASKAP_specline_val.py Eridanus_processed 1234
#
# Author: Bi-Qing For
# Email: biqing.for@icrar.org
# 
# Modified Date: 27 August 2019
#
################################################################################

import os.path
import sys
import glob
from astropy.io import fits
import numpy as np
from datetime import datetime
import PIL 
from PIL import Image
from astropy.stats import median_absolute_deviation
from astropy.utils.exceptions import AstropyWarning
import warnings
import math
from scipy.constants import k as k_B
import subprocess 
import matplotlib.pyplot as plt 
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
from astropy.io.fits import getheader

PI = math.pi

#ignore annoying astropy warnings 
warnings.simplefilter('ignore', category=AstropyWarning)

fig_dir = 'Figures'
main_dir = sys.argv[1]
sbid = str(sys.argv[2])

#################################################
#       Functions for the main program
#################################################

def get_FitsHeader(fitsimage):
    """
    Getting basic information of processed mosaic
    """
    hdu = getheader(fitsimage)
    bmaj = round(float(hdu['BMAJ'])*3600., 2)  #arcsec
    bmin = round(float(hdu['BMIN'])*3600., 2)  #arcsec

    return bmaj, bmin

def get_HIPASS (ra, dec):
    """
    Getting the HIPASS sources (HICAT; Meyer et al. 2004) within the 6x6 sq deg field through VizieR. 
    """

    print ("Retrieving HIPASS sources from Vizier. Depending on connection, this might take a while....")
    
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=['HIPASS', '_RAJ2000', '_DEJ2000', 'RVsp', 'Speak', 'Sint', 'RMS', 'Qual'], catalog = 'VIII/73/hicat')

    TOKS_RA = ra.split(":")
    ra_hr = float(TOKS_RA[0])
    ra_min = float(TOKS_RA[1])
    ra_sec = float(TOKS_RA[2])
    ra_deg = round(15.0*(ra_hr + ra_min/60. + ra_sec/3600.), 5) #Converting it to decimal degree
    TOKS_DEC = dec.split(".", 2)
    dec_deg = float(TOKS_DEC[0])
    dec_arcmin = float(TOKS_DEC[1])
    dec_arcsec = float(TOKS_DEC[2])
    dec_tdeg = round(dec_deg + dec_arcmin/60. + dec_arcsec/3600., 5) #Converting it to decimal degree
    
    hipass_result = v.query_region(coord.SkyCoord(ra=ra_deg, dec=dec_tdeg, unit=(u.deg, u.deg), frame='icrs'), width=[6*u.deg])

    hipass_cat = 'hipass.txt'
    print (hipass_result['VIII/73/hicat'], file=open(fig_dir + '/' + hipass_cat,'w'))

    return hipass_cat

    
def get_Version (param):
    """
    Getting the latest ASKAPsoft version used for the data reduction.
    """

    line = subprocess.check_output(['tail', '-5', param]) # Grab the last 5 lines
    str_line = line.decode('utf-8')
    newline = str_line.splitlines()[0] # This picks up the first line of the 5 
    TOKS = newline.split()
    askapsoft = TOKS[-1]

    return askapsoft
        

def get_Flagging (flagging_file):
    """
    Getting flagging statistics.
    """

    line = subprocess.check_output(['tail', '-1', flagging_file]) #Grab the last line
    str_line = line.decode('utf-8')
    TOKS = str_line.split()
    flagged_stat = TOKS[-1]

    return flagged_stat

    
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


def make_Thumbnail(image, thumb_img, sizeX, sizeY, dir):
    """
    Making thumbnail image.
    """
    
    size = sizeX, sizeY
    im = Image.open(image)
    im.thumbnail(size)
    im.save(dir+ '/' + thumb_img)

    return thumb_img


def cal_beam_MADMFD(infile):
    """
    Calculating the MAD of max flux density of each beam.
    """

    data = np.loadtxt(infile)
    maxfdensity = data[:,8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 2)

    return mad_maxfdensity
    

def cal_beam_AvgRMS(infile):
    """
    Calculating the average RMS of each beam. 
    """
    
    data = np.loadtxt(infile)
    rms = data[:,3]
    avg_rms = round(np.mean(rms), 2)
    
    return avg_rms
    

def cal_mosaic_Stats (infile):
    """
    Calculating MAD RMS and mean MADFM for the mosaic cube. 
    """

    data = np.loadtxt(infile)
    rms, madfm = data[:,3], data[:,5]
    
    med_madfm = np.median(madfm)
    mad_rms = round(median_absolute_deviation(rms), 2)

    return mad_rms, med_madfm

def qc_Bad_Chans (infile, med_madfm):
    """
    Checking for bad channels in the mosaic cube.
    """

    BAD_CHAN = []

    stat_file = open(infile, 'r')
    LINES = stat_file.readlines()[2:]
    stat_file.close()
    value = med_madfm + 0.4  # Deviation from the med_madfm. Need to check with larger sample of data to decide the best value. 

    for i in range(len(LINES)):
        line = LINES[i]
        TOKS = line.split()
        chan = TOKS[0]
        madfm = float(TOKS[5])
        if madfm > value:
            BAD_CHAN.append(chan)

    if BAD_CHAN == []:
        BAD_CHAN.append('none')
        QC_badchan_id = 'good'
    else:
        QC_badchan_id = 'bad'
        
    return BAD_CHAN, QC_badchan_id
    

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

def cal_binnedAvg(dataArray, N):
    """
    Calculating the average value of each bin. N is the specified bin number.  
    """
    
    mean_bin = np.cumsum(dataArray, 0)[N-1::N]/float(N)
    mean_bin[1:] = mean_bin[1:] - mean_bin[:-1]
    return mean_bin

"""
def qc_Max_Flux_Density (infile, delta_freq_range, mean_beamMADMFD):

    Evaluating the max flux density (MFD) of the mosaic.
    Metric to check for the effect of solar interference. 
    First, check the processed mosaic contains at least 5 MHz of data. 
    Then, calculating the MAD of MFD and comparing it to the mean of MADMFD of central 16 beams. 


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

    
def qc_RMS(infile, theoretical_rms_mjy):

    Evaluationg the RMS values of mosaic in ~1 MHz interval.  

    
    N = 54
    count = 0
    
    data = np.loadtxt(infile)
    rms = data[:,3]
    bin_rms = cal_binnedAvg(rms, N)
    if bin_rms.any > theoretical_rms_mjy:
        count += 1  
    
    return bin_rms

"""

def BeamStat_plot(infile, item):
    """
    Plotting and visualising statistics of 36 beams. 
    """

    if item == 'MADMFD':
        vmin = 0.1   # This is a conservative cut off based on M83 field. 0.5 mJy/beam is more realistic.
        vmax = 1.0
        title = 'MAD Max Flux Density'
        plot_name = 'beamStat_MADMFD.png'
        saved_fig = fig_dir+'/'+plot_name
        

    if item == 'Avg_RMS':
        vmin = 2.0
        vmax = 4.0
        title = 'Mean RMS'
        plot_name = 'beamStat_AvgRMS.png'
        saved_fig = fig_dir+'/'+plot_name
    
    n = [26,25,24,23,22,21,27,10,9,8,7,20,28,11,3,1,6,19,29,12,2,0,5,18,30,13,14,15,4,17,31,32,33,34,35,16] # beam number

    XPOS, YPOS = [], []

    x=0
    for j in range(0,6,1):
        x += 0.1
        y=0
        for k in range(0,6,1):
            y += 0.2
            XPOS.append(x)
            YPOS.append(y)

    for i in range(36):
        bnum = n[i]
        if bnum < 10:
            infile = 'cubeStats-image.restored.i.SB'+ sbid +'.cube.'+ field + '.beam0' + str(bnum) +'.contsub.txt'
        else:
            infile = 'cubeStats-image.restored.i.SB'+ sbid +'.cube.'+ field + '.beam' + str(bnum) +'.contsub.txt'
        if os.path.isfile(main_dir + '/' + field +'/'+infile):
            if item == 'MADMFD': 
                beamstat = cal_beam_MADMFD(main_dir + '/' + field +'/'+infile)
            if item == 'Avg_RMS':
                beamstat = cal_beam_AvgRMS(main_dir + '/' + field +'/'+infile)
#        else:
#            beamstat = 0.5

        plt.scatter([XPOS[i]], [YPOS[i]], s=1500, c=[beamstat], cmap='RdYlGn_r', edgecolors='black', vmin=vmin, vmax=vmax)
        plt.text(XPOS[i], YPOS[i], n[i])

    plt.xlim(0,0.7)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    cb.set_label('mJy / beam')
    plt.savefig(saved_fig)
    plt.close()

    return saved_fig, plot_name



"""
def check_Cleaning ():
    
    Checking if cleaning has performed down to the defined threshold. Both major and minor cycles. 
    
    This will involve ncore specified. 
"""    
"""
def plot(infile, x, y, c=None, yerr=None, figure=None, arrows=None, xlabel='', ylabel='')
"""


###########################################################
# Main program where it calls all the functions
###########################################################

# Required files 

metafile = sorted(glob.glob('metadata/mslist-*txt'))[0]
metafile_science = sorted(glob.glob('metadata/mslist-scienceData*txt'))[0]
cubestat_contsub = glob.glob(main_dir + '/cubeStats*contsub.txt')[0]
flagging_file = glob.glob('slurmOutput/flag.out.txt')[0]  # temporary file from Wasim
param_file = sorted(glob.glob('slurmOutput/*.sh'))
fitsimage = (main_dir+'/image.restored.i.SB' + sbid + '.cube.contsub.fits')


# Check if there is more than one parameter input .sh file in the slurmOutput directory.
# If it does, select the latest one.
# More than one version is used. Reporting the latest version of ASKAPSoft for final data reduction. 

if len(param_file) >=1:
    index = len(param_file)
    param = param_file[index-1]
else:
    param = param_file[0]


#############################    
# HTML related tasks
#############################

html_name = 'spectral_report_SB' + sbid + '.html'

### Making thumbnail images

sizeX = 70
sizeY = 70

cube_plots = glob.glob(main_dir + '/cubePlot*.png')
beamNoise_plots = glob.glob(main_dir + '/beamNoise*.png')
beamMinMax_plots = glob.glob(main_dir + '/beamMinMax*.png')

thumb_cubeplots = []
thumb_beamNoise = []
thumb_beamMinMax = []

for image in cube_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.partition('/')[2]
    thumb_cubeplots.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, main_dir)

for image in beamNoise_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.partition('/')[2]
    thumb_beamNoise.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, main_dir)

for image in beamMinMax_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.partition('/')[2]
    thumb_beamMinMax.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, main_dir)


# Calling the functions

askapsoft = get_Version(param)
n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec = get_Metadata(metafile)
chan_width, bw, cfreq = get_Metadata_freq(metafile_science)
start_freq, end_freq = get_Frequency_Range(cubestat_contsub)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.
hipass_cat = get_HIPASS(ra, dec)
bmaj, bmin = get_FitsHeader(fitsimage)

tobs_hr = round(tobs/3600.,2) # convert tobs from second to hr
chan_width_kHz = round(chan_width/1000.,3) # convert Hz to kHz

freq_range = str(start_freq)+'--'+str(end_freq)
delta_freq_range = end_freq - start_freq


### Validation parameters

flagged_stat = get_Flagging(flagging_file)

beamstat_contsub = glob.glob(main_dir + '/' + field +'/cubeStats*beam??.contsub.txt')

if not os.path.isdir(fig_dir):
    os.system('mkdir '+ fig_dir)

beam_MADMFD_fig, MADMFD_plot = BeamStat_plot(cubestat_contsub, 'MADMFD')
thumb_img = 'thumb_'+ str(sizeX) + '_'+ MADMFD_plot
make_Thumbnail(beam_MADMFD_fig, thumb_img, sizeX, sizeY, fig_dir)

beam_Avg_RMS_fig, AvgRMS_plot = BeamStat_plot(cubestat_contsub, 'Avg_RMS')
thumb_img = 'thumb_'+ str(sizeX) + '_'+ AvgRMS_plot
make_Thumbnail(beam_Avg_RMS_fig, thumb_img, sizeX, sizeY, fig_dir)


#QC_mad_maxfden, QC_maxfden_id = qc_Max_Flux_Density(cubestat_contsub, delta_freq_range) #Continuum subtracted
QC_mad_maxfden = '2.0'
QC_maxfden_id = 'good'

mad_rms, med_madfm = cal_mosaic_Stats(cubestat_contsub)
bad_chans, QC_badchan_id = qc_Bad_Chans(cubestat_contsub, med_madfm)

#bin_value = qc_RMS(cubestat_contsub, theoretical_rms_mjy)



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
                        <th>ASKAPsoft<br>version*</th>
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
                    <th>Total <br> Flagged Visibilities</th>
                    <th>MAD RMS <br> mJy/beam</th>
                    <th>Expected RMS <br> mJy/beam</th>
                    <th>Bad Channel</th>
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
                    <td>{12}</td>
                    <td>{13}</td>
                    <td>{14:.1f}</td>
                    <td id='{15}'>{16}
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
                               flagged_stat,
                               mad_rms,
                               theoretical_rms_mjy,
                               QC_badchan_id,
                               bad_chans))


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


html.write("""</td>
                    </tr>
                    </table>
                    <table width="100%" border="1">
                    <tr>
                    <th colspan="2">Continuum Subtracted Beam Cube</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>MAD Max Flux Density</br></p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>Mean RMS</br></p>
                    """.format(fig_dir+'/'+'beamStat_MADMFD.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_MADMFD.png',
                               sizeX,
                               sizeY,
                               fig_dir+'/'+'beamStat_AvgRMS.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_AvgRMS.png',
                               sizeX,
                               sizeY))

html.write("""</td> 
              </tr>
              </table>
              <h2 align="middle">Validation Metrics</h2>
              <table width="100%" border="1">
              <tr>
              <th>MAD Max Flux Density <br>(mJy/beam)</th>
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
              <h2 align="middle">HIPASS sources within 6x6 sq degree</h2>
              <table>
              <tr align="middle">
              <td id="cell">
              <form>
              <input type="button" value="Click here" onclick="window.location.href='{0}'" style="font-size:20px; width=50%; height=50%"</>
              </form>
              """.format(fig_dir+'/'+hipass_cat))

### Finish HTML report with generated time stamp

html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>* If more than one version of ASKAPsoft is used for the whole reduction, the latest one is reported.<br>
              Generated at {0} <br>
              <i> Report bugs to 
              <a href="mailto:biqing.for@icrar.org">Bi-Qing For</a></i>
              </p>

              </body>
              </html>""".format(str(datetime.now())))





html.close()
print ("Spectral line validation report written to '{0}'.".format(html_name))
                    
