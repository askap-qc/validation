################################################################################
# This script generates an ASKAP spectral line cube HTML report. 
#
# Compatibility: Python version 3.x
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- metadata/mslist-*101.txt
# -- slurmOutput/*sh
# -- image.restored.i.SB<SBID>.cube.contsub.fits
# -- diagnostics/cubestats-<field>/*txt
# -- diagnostics/*png
#
# Output files: all saved in Figures directory.
#
# To run type:
#
# "python <script name> <SBID>
# e.g. python ASKAP_specline_val.py 8170
#
# Author: Bi-Qing For
# Email: biqing.for@icrar.org
# 
# Modified Date: 4 Sept 2019
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
import warnings
from astropy.utils.exceptions import AstropyWarning
from scipy.optimize import OptimizeWarning
import math
from scipy.constants import k as k_B
import subprocess 
import matplotlib.pyplot as plt 
import astropy.coordinates as coord
import astropy.units as u
from astropy.io.fits import getheader
from scipy.stats import iqr
import matplotlib.pylab as pylab
from scipy.optimize import curve_fit
from numpy import inf
from scipy import asarray as ar,exp
import matplotlib.patches as mpatches

PI = math.pi

#################################################
#       Functions for the main program
#################################################

def BeamPosition():
    """
    Defining 36 beam positions for plotting.
    """
    
    XPOS, YPOS = [], []

    x=0
    for j in range(0,6,1):
        x += 0.1
        y=0
        for k in range(0,6,1):
            y += 0.2
            XPOS.append(x)
            YPOS.append(y)

    return XPOS, YPOS

def gauss(x, *p):
    """ 
    Fitting a Gaussian function.
    """
    A, mu, sigma = p

    return A*exp(-(x-mu)**2/(2.*sigma**2))


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
    from astroquery.vizier import Vizier

    print ("Retrieving HIPASS sources from Vizier. Depending on server connection, this might take a while....")
    
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

    nBlocks = 6  # these are the number of correlator cards (PILOT survey value)
    
    obs_date = 'Observed from'
    code = 'Code'
    duration = 'Total elapsed time'
    antenna = 'antennas'
    frame = 'Frame'
    
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
        if line.find(frame) >= 0:
            next_line = LINES[i+1]
            TOKS = next_line.split()
            total_obs_bw = float(TOKS[10])*nBlocks/1000.0 # kHz to MHz 
            
    return n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, total_obs_bw

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
            cfreq = TOKS[12] #MHz
            nchan = TOKS[7]

    return chan_width, cfreq, nchan    
    

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

    mosaic_bad_chan = 'mosaic_badchans.txt'
    print (','.join(BAD_CHAN), file=open(fig_dir + '/' + mosaic_bad_chan,'w'))

    n_bad_chan = len(BAD_CHAN)

    # Check if number of bad channel recorded is 1. If yes, check if is it a none keyword.
    # If yes, number of bad channel should be 0.
    
    if n_bad_chan == 1:
        with open(fig_dir + '/' + mosaic_bad_chan) as f:
            if 'none' in f.read():
                n_bad_chan = 0
                print ('yes')
    
    return n_bad_chan, mosaic_bad_chan, QC_badchan_id
    

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


def qc_NoiseRank(spread):
    """
    Evaluating the 1-percentile noise rank distribution. Sigma of the Gaussian function is used as a metric. 
    """

    sigma = 0.5 # for a Gaussian
    
    if (spread <= 2.5*sigma):
        qc_label = 'good'
    elif (2.5*sigma < spread < 4*sigma):
        qc_label = 'ok'
    else:
        qc_label = 'bad'

    return qc_label

def NoiseRank_histplot(nchan):
    
    ID_LABEL = []
    plot_name = 'beam_1pctile_hist_SB'+ sbid + '.png'
    saved_fig = fig_dir + '/' + plot_name
    file_dir = diagnostics_dir +'/cubestats-'+ field 
    basename = '/cubeStats-image.restored.i.SB'+ sbid +'.cube.'+ field

    params = {'axes.labelsize': 6,
              'axes.titlesize':6,
              'xtick.labelsize':5,
              'ytick.labelsize':5,
              'font.size':6}

    pylab.rcParams.update(params)

    fig, axs = plt.subplots(6,6)
    fig.subplots_adjust(hspace = .004, wspace=.004, bottom=0.05)
    fig.text(0.5, 0.007, '1-percentile noise level (mJy / beam)', ha='center')
    fig.text(0.01, 0.5, 'log N', ha='center', rotation='vertical')

    axs = axs.ravel()

    for i in range(0,36):
        infile = file_dir + basename +'.beam%02d.contsub.txt'%(i)
        data = np.loadtxt(infile)
        onepctile = data[:,6]
        median_val = np.median(onepctile)
        # if statement is needed to rule out really bad data without having to do the Gaussian fitting
        if (median_val > 1000.0) or (median_val < -1000.0):
            ID_LABEL.append('bad')
            axs[i].set_xlim(-1, 1)
            axs[i].set_ylim(0, 3)
            axs[i].title.set_text('Beam%02d' %(i))
        else:
            upper_range = median_val + 3
            lower_range = median_val - 3
            x = onepctile[(onepctile < upper_range) & (onepctile > lower_range)]  # exclude outliers
            xmax_val = np.max(x) 
            xmin_val = np.min(x)

            # Freedman-Diaconis rule. Nchan includes all processed channels, not excluding outliers. 
            bin_width = 2*iqr(x)*nchan**(-1/3) 
            n_bins = int((xmax_val - xmin_val)/bin_width)
    
            hist, bins = np.histogram(onepctile, bins=n_bins, range=(xmin_val-3, xmax_val+3))
            with np.errstate(divide='ignore'):  # ignore division of zero 
                N = np.log10(hist)   # get log N for y-axis
                N[N == -inf] = 0

            xcenter = (bins[:-1] + bins[1:]) / 2
            ymax_val = np.max(N)
            median_val_x = np.median(x)
            
            # Fitting a Gaussian and use spread (sigma) as a metric
            guess=[ymax_val, median_val_x, 5.0]
            coeff, var_matrix = curve_fit(gauss, xcenter, N, guess)
            spread = round(np.abs(coeff[2]), 1)
            ID_LABEL.append(qc_NoiseRank(spread))
    
            axs[i].bar(xcenter, N)
            axs[i].plot(xcenter,gauss(xcenter,*coeff),'r-',lw=1)    
            axs[i].set_xlim(xmin_val-3, xmax_val+3)
            axs[i].set_ylim(0, ymax_val+3)
            axs[i].title.set_text('Beam%02d' %(i))

    plt.tight_layout()
    plt.savefig(saved_fig)
    plt.close()

    return saved_fig, plot_name, ID_LABEL


def NoiseRank_QCplot(list_id_label, n):


    legend_dict = { 'good' : '#00FF00', 'ok' : '#FFD700', 'bad' : '#CD5C5C' }
    patchList = []
    XPOS, YPOS = BeamPosition()

    plot_name = 'beam_1pctile_qc_SB' + sbid + '.png'
    saved_fig = fig_dir + '/' + plot_name
    
    for i in range(36):
        bnum = n[i]
        if (list_id_label[bnum] =='good'):
            color_code = '#00FF00'
        elif (list_id_label[bnum] =='ok'): 
            color_code = '#FFD700'
        else:
            color_code = '#CD5C5C'
        
        plt.scatter([XPOS[i]], [YPOS[i]], s=1500, color=color_code, edgecolors='black')
        plt.text(XPOS[i], YPOS[i], n[i])

    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    plt.legend(handles=patchList)
    plt.xlim(0,0.8)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title('1-percentile noise rank')
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


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

def BeamStat_plot(item, n):
    """
    Plotting and visualising statistics of 36 beams. 
    """
    file_dir = diagnostics_dir +'/cubestats-'+ field 
    basename = '/cubeStats-image.restored.i.SB'+ sbid +'.cube.'+ field  
    
    if item == 'MADMFD':
        vmin = 0.5   # This is a conservative cut off based on M83 field. 0.5 mJy/beam is more realistic.
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
    
    beamXPOS, beamYPOS = BeamPosition()
    
    for i in range(36):
        bnum = n[i]
        infile = file_dir + basename +'.beam%02d.contsub.txt'%(bnum)
        if os.path.isfile(infile):
            if item == 'MADMFD': 
                beamstat = cal_beam_MADMFD(infile)
            if item == 'Avg_RMS':
                beamstat = cal_beam_AvgRMS(infile)
#        else:
#            beamstat = 0.5

        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1500, c=[beamstat], cmap='RdYlGn_r', edgecolors='black', vmin=vmin, vmax=vmax)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

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

#ignore astropy warnings 
warnings.simplefilter('ignore', AstropyWarning)   

fig_dir = 'Figures'
sbid = str(sys.argv[1])
n = [26,25,24,23,22,21,27,10,9,8,7,20,28,11,3,1,6,19,29,12,2,0,5,18,30,13,14,15,4,17,31,32,33,34,35,16] # beam number
html_name = 'spectral_report_SB' + sbid + '.html'
diagnostics_dir = 'diagnostics'

if not os.path.isdir(fig_dir):
    os.system('mkdir '+ fig_dir)


# Required files 

metafile = sorted(glob.glob('metadata/mslist-*txt'))[0]
metafile_science = sorted(glob.glob('metadata/mslist-scienceData*txt'))[0]
flagging_file = glob.glob('slurmOutput/flag.out.txt')[0]  # temporary file from Wasim
param_file = sorted(glob.glob('slurmOutput/*.sh'))
fitsimage = ('image.restored.i.SB' + sbid + '.cube.contsub.fits')


# Check if there is more than one parameter input .sh file in the slurmOutput directory.
# If it does, select the latest one.
# If more than one version is used. Reporting the latest version of ASKAPSoft for final data reduction. 

if len(param_file) >=1:
    index = len(param_file)
    param = param_file[index-1]
else:
    param = param_file[0]


n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, total_obs_bw = get_Metadata(metafile)
askapsoft = get_Version(param)
chan_width, cfreq, nchan = get_Metadata_freq(metafile_science)
tobs_hr = round(tobs/3600.,2) # convert tobs from second to hr
chan_width_kHz = round(chan_width/1000.,3) # convert Hz to kHz

cubestat_linmos_contsub = glob.glob(diagnostics_dir+ '/cubestats-' + field + '/cubeStats*linmos.contsub.txt')[0] #mosaic contsub statistic

start_freq, end_freq = get_Frequency_Range(cubestat_linmos_contsub)
freq_range = str(start_freq)+'--'+str(end_freq)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.
bmaj, bmin = get_FitsHeader(fitsimage)
delta_freq_range = end_freq - start_freq
flagged_stat = get_Flagging(flagging_file)
mad_rms, med_madfm = cal_mosaic_Stats(cubestat_linmos_contsub)
n_bad_chan, mosaic_bad_chan, QC_badchan_id = qc_Bad_Chans(cubestat_linmos_contsub, med_madfm)
hipass_cat = get_HIPASS(ra, dec)

    
#############################    
# HTML related tasks
#############################

### Making thumbnail images

sizeX = 70
sizeY = 70

cube_plots = glob.glob(diagnostics_dir + '/cubestats-' + field + '/*linmos*.png')  #Mosaic statistic
beamNoise_plots = glob.glob(diagnostics_dir + '/beamNoise*.png') #beam-by-beam statistic
beamMinMax_plots = glob.glob(diagnostics_dir +'/beamMinMax*.png') #beam-by-beam statistic

thumb_cubeplots = []
thumb_beamNoise = []
thumb_beamMinMax = []

for image in cube_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.split('/')[2]
    thumb_cubeplots.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

for image in beamNoise_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.split('/')[1]
    thumb_beamNoise.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

for image in beamMinMax_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.split('/')[1]
    thumb_beamMinMax.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)


beam_MADMFD_fig, MADMFD_plot = BeamStat_plot('MADMFD', n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ MADMFD_plot
make_Thumbnail(beam_MADMFD_fig, thumb_img, sizeX, sizeY, fig_dir)

beam_Avg_RMS_fig, AvgRMS_plot = BeamStat_plot('Avg_RMS', n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ AvgRMS_plot
make_Thumbnail(beam_Avg_RMS_fig, thumb_img, sizeX, sizeY, fig_dir)

#QC_mad_maxfden, QC_maxfden_id = qc_Max_Flux_Density(cubestat_contsub, delta_freq_range) #Continuum subtracted
#QC_mad_maxfden = '2.0'
#QC_maxfden_id = 'good'

#bin_value = qc_RMS(cubestat_contsub, theoretical_rms_mjy)

beam_1pctile_histfig, onepctile_plot, list_id_label = NoiseRank_histplot(float(nchan))
thumb_img = 'thumb_' + str(sizeX) + '_' + onepctile_plot
make_Thumbnail(beam_1pctile_histfig, thumb_img, sizeX, sizeY, fig_dir)

beam_1pctile_QCfig, onepctile_QCplot = NoiseRank_QCplot(list_id_label, n)
thumb_img = 'thumb_' + str(sizeX) + '_' + onepctile_QCplot
make_Thumbnail(beam_1pctile_QCfig, thumb_img, sizeX, sizeY, fig_dir)

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

                  #ok {
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
                        <th>Total Bandwidth <br>(MHz)</th>
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
                        <td>{8}""".format(sbid,
                                          n_ant,
                                          start_obs_date,
                                          end_obs_date,
                                          tobs_hr,
                                          field,
                                          ra,
                                          dec,
                                          total_obs_bw))
                                          

html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Processed Image Cube</h2>
                    <table width="100%" border="1">
                    <tr>
                        <th>ASKAPsoft<br>version*</th>
                        <th>Frequency Range<br>(MHz)</th>
                        <th>Central Frequency<br>(MHz)</th>
                        <th>Channel Width<br>(kHz)</th>
                        <th>Synthesised Beam bmaj<br>(arcsec)</th>
                        <th>Synthesised Beam bmin<br>(arcsec)</th>
                    </tr>
                    <tr align="middle">
                        <td>{0}</td>
                        <td>{1}</td>
                        <td>{2}</td>
                        <td>{3}</td>
                        <td>{4}</td>
                        <td>{5}""".format(askapsoft,
                                          freq_range,
                                          cfreq,
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
                    <th>Number of Bad Channel</th>
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
                    <form action="{17}" method="get" target="_blank">
                     <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
                    </form>
                    """.format(cube_plots[0],
                               fig_dir+'/' + thumb_cubeplots[0],
                               sizeX,
                               sizeY,
                               cube_plots[1],
                               fig_dir+'/'+ thumb_cubeplots[1],
                               sizeX,
                               sizeY,
                               cube_plots[2],
                               fig_dir+'/'+ thumb_cubeplots[2],
                               sizeX,
                               sizeY,
                               flagged_stat,
                               mad_rms,
                               theoretical_rms_mjy,
                               QC_badchan_id,
                               n_bad_chan,
                               fig_dir+'/' + mosaic_bad_chan))


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
                    <br><p>Min, Max, 1 percentile</p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile</p>
                    </td>
                    <td>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile</p>
                    </td>
                    <tr align="middle">
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</p>
                    </td>
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</p>
                    </td>
                    <td>
                    <a href="{20}" target="_blank"><img src="{21}" width="{22}" height="{23}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</p>
                    """.format(beamMinMax_plots[1],
                               fig_dir+'/'+ thumb_beamMinMax[1],
                               sizeX,
                               sizeY,
                               beamMinMax_plots[0],
                               fig_dir+'/'+ thumb_beamMinMax[0],
                               sizeX,
                               sizeY,
                               beamMinMax_plots[2],
                               fig_dir+'/'+ thumb_beamMinMax[2],
                               sizeX,
                               sizeY,
                               beamNoise_plots[0],
                               fig_dir+'/'+ thumb_beamNoise[0],
                               sizeX,
                               sizeY,
                               beamNoise_plots[1],
                               fig_dir+'/'+ thumb_beamNoise[1],
                               sizeX,
                               sizeY,
                               beamNoise_plots[2],
                               fig_dir+'/'+ thumb_beamNoise[2],
                               sizeX,
                               sizeY))


html.write("""</td>
                    </tr>
                    </table>
                    <table width="100%" border="1">
                    <tr>
                    <th colspan="3">Continuum Subtracted Beam Cube</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>MAD Max Flux Density</p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>Mean RMS</p>
                    </td>
                    <td>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p>1-percentile noise rank</p>
                    """.format(fig_dir+'/'+'beamStat_MADMFD.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_MADMFD.png',
                               sizeX,
                               sizeY,
                               fig_dir+'/'+'beamStat_AvgRMS.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_AvgRMS.png',
                               sizeX,
                               sizeY,
                               beam_1pctile_histfig,
                               fig_dir+'/'+ 'thumb_' + str(sizeX) + '_' + onepctile_plot,
                               sizeX,
                               sizeY,
                               beam_1pctile_QCfig, 
                               fig_dir+'/'+ 'thumb_' + str(sizeX) + '_' + onepctile_QCplot,
                               sizeX,
                               sizeY))

###html.write("""</td> 
###              </tr>
###              </table>
###              <h2 align="middle">Validation Metrics</h2>
###              <table width="100%" border="1">
###              <tr>
###              <th>MAD Max Flux Density <br>(mJy/beam)</th>
###              <th>item</th>
###              <th>item</th>   
###              </tr>
###              <tr align="middle">
###              <td id='{0}'>{1}</td>
###              <td id='{2}'>{3}</td>
###              <td id='{4}'>{5}""".format(QC_maxfden_id,
###                                         QC_mad_maxfden,
###                                         QC_maxfden_id,
###                                         QC_mad_maxfden,
###                                         QC_maxfden_id,
###                                         QC_mad_maxfden))

html.write("""</td> 
              </tr>
              </table>
              <h2 align="middle">HIPASS sources within 6x6 sq degree**</h2>
              <table>
              <tr align="middle">
              <td id="cell">
              <form action="{0}" method="get" target="_blank">
                 <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
              </form>
              """.format(fig_dir+'/' + hipass_cat))

### Finish HTML report with generated time stamp

html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              * If more than one version of ASKAPsoft is used for the whole reduction, the latest one is reported.<br>
              ** Does not take into account field rotation. <br>  
              Generated at {0} <br>
              <i> Report bugs to 
              <a href="mailto:biqing.for@icrar.org">Bi-Qing For</a></i>
              </p>

              </body>
              </html>""".format(str(datetime.now())))

html.close()
print ("Spectral line validation report written to '{0}'.".format(html_name))
                    
