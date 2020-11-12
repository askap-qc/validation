################################################################################
#
# ASKAP-FLASH data validation script
#
# This script generates an ASKAP-FLASH spectral line cube HTML report. 
# The original script is written for ASKAP-WALLABY spectral line cube by Bi-Qing For.
#
# Compatibility: Python version 3.x
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- metadata/mslist-cal*txt
# -- metadata/mslist-*101.txt
# -- slurmOutput/*sh
# -- image.restored.i.SB<SBID>.cube.contsub.fits (optional see comment)
# -- diagnostics/cubestats-<field>/*txt
# -- diagnostics/*png
# -- diagnostics/Flagging_Summaries/*SL.ms.flagSummary 
# -- SpectralCube_BeamLogs/*.txt
#
# Output files: All saved in Figures directory. 
#
# To run type:
#
# python <script name> <SBID>
# e.g. python ASKAP_specline_val.py 13293
#
# Modified for FLASH: Hyein Yoon
# Email: hyein.yoon [at] sydney.edu.au
# Last modified: 12 November 2020
#
# Original author for WALLABY: Bi-Qing For
# Email: biqing.for [at] icrar.org
# Modified Date: 24 March 2020
#
################################################################################


################################################################################
# Modules for the script
################################################################################

import os
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
import astropy.coordinates as coord
import astropy.units as u
from astropy.io.fits import getheader
from scipy.stats import iqr
from scipy.optimize import curve_fit
from numpy import inf
from scipy import asarray as ar,exp
from astroquery.vizier import Vizier


# This step is necessary to avoid matplotlib using the Xwindows backend. 
# Otherwise, it does not work on galaxy. 

import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('No display is found. Using the non-interactive Agg backend.')
    mpl.use('Agg')
import matplotlib.pyplot as plt 
#        else:
#            beamstat = 0.5
import matplotlib.pylab as pylab
import matplotlib.patches as mpatches

PI = math.pi


################################################################################
# Functions for the main program
################################################################################

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


################################################################################
# Getting the FIRST sources (2014Dec17; Helfand+ 2015) within the 6x6 sq deg field through VizieR
################################################################################

def get_FIRST(ra, dec):
    from astroquery.vizier import Vizier # for installation: conda install -c astropy astroquery

    print ("Retrieving FIRST sources from Vizier. Depending on server connection, this might take a while....")
    
#   Vizier.ROW_LIMIT = 1000 # Note: it doesn't work!
#   v = Vizier(columns=['FIRST', '_RAJ2000', '_DEJ2000', 'p(S)', 'Fpeak', 'Fint', 'Maj', 'Min'], catalog = 'VIII/92/first14')
    v = Vizier(columns=['all'], catalog = 'VIII/92/first14')
    Vizier.ROW_LIMIT = -1

    TOKS_RA = ra.split(":")
    ra_hr = float(TOKS_RA[0])
    ra_min = float(TOKS_RA[1])
    ra_sec = float(TOKS_RA[2])
    ra_deg = round(15.0*(ra_hr + ra_min/60. + ra_sec/3600.), 5) # Converting it to decimal degree
    TOKS_DEC = dec.split(".", 2)
    dec_deg = float(TOKS_DEC[0])
    dec_arcmin = float(TOKS_DEC[1])
    dec_arcsec = float(TOKS_DEC[2])
    dec_tdeg = round(dec_deg + dec_arcmin/60. + dec_arcsec/3600., 5) # Converting it to decimal degree
    
    first_result = v.query_region(coord.SkyCoord(ra=ra_deg, dec=dec_tdeg, unit=(u.deg, u.deg), frame='icrs'), width=[6*u.deg])

    print(first_result)
    first_cat = 'first.txt'
    print (first_result['VIII/92/first14'],file=open(fig_dir + '/' + first_cat,'w'))

    return first_cat


################################################################################
# Getting the NVSS sources (Condon+ 1998) within the 6x6 sq deg field through VizieR. 
################################################################################
  
def get_NVSS(ra, dec):
    from astroquery.vizier import Vizier # for installation: conda install -c astropy astroquery
 
    print ("Retrieving NVSS sources from Vizier. Depending on server connection, this might take a while....")
    
    v = Vizier(columns=['NVSS', '_RAJ2000', '_DEJ2000', 'S1.4', 'MajAxis', 'MinAxis'], catalog = 'VIII/65/nvss')
    v.ROW_LIMIT = -1

    TOKS_RA = ra.split(":")
    ra_hr = float(TOKS_RA[0])
    ra_min = float(TOKS_RA[1])
    ra_sec = float(TOKS_RA[2])
    ra_deg = round(15.0*(ra_hr + ra_min/60. + ra_sec/3600.), 5) # Converting it to decimal degree
    TOKS_DEC = dec.split(".", 2)
    dec_deg = float(TOKS_DEC[0])
    dec_arcmin = float(TOKS_DEC[1])
    dec_arcsec = float(TOKS_DEC[2])
    dec_tdeg = round(dec_deg + dec_arcmin/60. + dec_arcsec/3600., 5) # Converting it to decimal degree
    
    nvss_result = v.query_region(coord.SkyCoord(ra=ra_deg, dec=dec_tdeg, unit=(u.deg, u.deg), frame='icrs'), width=[6*u.deg])

    print(nvss_result) 
    nvss_cat = 'nvss.txt'
    print (nvss_result['VIII/65/nvss'], file=open(fig_dir + '/' + nvss_cat,'w'))

    return nvss_cat


"""
################################################################################  
# Getting the HIPASS sources (HICAT; Meyer et al. 2004) within the 6x6 sq deg field through VizieR. 
################################################################################

def get_HIPASS(ra, dec):

    print ("Retrieving HIPASS sources from Vizier. Depending on server connection, this might take a while......")
    
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=['HIPASS', '_RAJ2000', '_DEJ2000', 'RVsp', 'Speak', 'Sint', 'RMS', 'Qual'], catalog = 'VIII/73/hicat')

    TOKS_RA = ra.split(":")
    ra_hr = float(TOKS_RA[0])
    ra_min = float(TOKS_RA[1])
    ra_sec = float(TOKS_RA[2])
    ra_deg = round(15.0*(ra_hr + ra_min/60. + ra_sec/3600.), 5) # Converting it to decimal degree
    TOKS_DEC = dec.split(".", 2)
    dec_deg = float(TOKS_DEC[0])
    dec_arcmin = float(TOKS_DEC[1])
    dec_arcsec = float(TOKS_DEC[2])
    dec_tdeg = round(dec_deg + dec_arcmin/60. + dec_arcsec/3600., 5) # Converting it to decimal degree
    
    hipass_result = v.query_region(coord.SkyCoord(ra=ra_deg, dec=dec_tdeg, unit=(u.deg, u.deg), frame='icrs'), width=[6*u.deg])

    hipass_cat = 'hipass.txt'
    print (hipass_result['VIII/73/hicat'], file=open(fig_dir + '/' + hipass_cat,'w'))

    return hipass_cat
"""

def get_Version(param):
    """
    Getting the latest ASKAPsoft version that is used for the data reduction.
    """

    line = subprocess.check_output(['tail', '-5', param]) # Grab the last 5 lines
    str_line = line.decode('utf-8')
    newline = str_line.splitlines()[0] # This picks up the first line of the 5 
    TOKS = newline.split()
    askapsoft = TOKS[-1]

    return askapsoft
        

def get_Flagging_KeyValues(flagging_file):
    """
    Getting Flagging Key Values. 
    """

    flag_infile = open(flagging_file, 'r')
    LINES = flag_infile.readlines()[:6]
    flag_infile.close()
    
    N_Rec = 'nRec'  # Total number of spectra feeds into the synthesis image. This is not always constant so grab the value beam-by-beam.
    N_Chan = 'nChan'  # Total number of channel

    # Search for keywords in the file
    
    for i in range(len(LINES)):
        line = LINES[i]
        if line.find(N_Rec) >=0:
            TOKS = line.split()
            n_Rec = float(TOKS[2])
        if line.find(N_Chan) >=0:
            TOKS = line.split()
            n_Chan = float(TOKS[2])

    exp_count = n_Rec*35 #counting antenna number from zero based on the recorded data
    
    return n_Rec, n_Chan, exp_count


def get_Flagging(flagging_file, n_Rec, nChan, exp_count):
    """
    Getting flagging statistics and finding out beam-by-beam antenna based (completely) flagging. 
    """

    line = subprocess.check_output(['tail', '-1', flagging_file]) #Grab the last line
    str_line = line.decode('utf-8')
    TOKS = str_line.split()
    total_flagged_pct = float(TOKS[-2]) #data+autocorrelation
    total_uv = float(TOKS[7])

    # Getting data flagged percentage
    
    autocorr_flagged_pct = (36 * n_Rec * n_Chan / total_uv)*100.0
    data_flagged_pct = round(total_flagged_pct - autocorr_flagged_pct, 3)

    # Finding out which antenna has been flagged completely.
    ANT1, ANT2, FLAG = [], [], [] 
    with open(flagging_file, 'r') as f:
        for line in f:
            if "#" not in line:  # grep -v "#"
                if "Flagged" not in line:   # grep -v "Flagged"
                    TOKS=line.split()       
                    ant1 = int(TOKS[3])
                    ant2 = int(TOKS[4])
                    flag = float(TOKS[6])
                    if (ant1 < ant2) and (flag == 100): # extract non-correlated antenna pairs with 100 percent flagging
                        ANT1.append(ant1)
                        ANT2.append(ant2)
                        FLAG.append(flag)

    with open('temp.dat' ,'w') as newlist:
        for j in range(len(ANT1)):
            newlist.write('%i %i %.2f\n' %(ANT1[j],ANT2[j],FLAG[j]))
        
    newlist.close()
    f.close()

    data = np.genfromtxt('temp.dat')
    ant1, ant2, flag = data[:,0], data[:,1], data[:,2]

    ANT_NAME = []
    for x in range(0,36):
        count1 = np.count_nonzero(ant1 == x)
        count2 = np.count_nonzero(ant2 == x)
        total_count = count1 + count2
        if total_count == exp_count:
            ant_num = x+1
            ant_name = 'ak'+ str(ant_num)
            ANT_NAME.append(ant_name)

    total_flagged_ant = len(ANT_NAME)
       
    os.system('rm temp.dat')

    flag_ant_file = 'flagged_antenna.txt'
    ffile = open(fig_dir + '/'+ flag_ant_file,'a')
    
    if total_flagged_ant > 1:
        ffile.write(flagging_file[-24:-18])
        ffile.write('\n')
        for item in ANT_NAME:
            ffile.write(item)
            ffile.write('\n')
    else:
        ffile.write(flagging_file[-24:-18])
        ffile.write('\n none \n')

    ffile.close()
                
    return data_flagged_pct, total_flagged_ant, flag_ant_file


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
    Frequency range and total channel of the mosaic
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

def qc_Missing_Data(infile):
    """
    Checking for missing data channels in the mosaic cube. 
    This helps checking for the correlator drop out or insufficient RFI flagging. 
    e.g. 4MHz (correlator block) ~ 216 channels. 
    """

    data = np.loadtxt(infile)
    maxfdensity = data[:,8]

    count = sum(map(lambda x : x==0, maxfdensity))

    # Not using 216 channels because Doppler correction can affect the number
    if count >= 100:
        QC_mdata_id = 'bad'
        QC_mdata_keyword = 'Yes > 100'
    elif (0 < count < 100):
        QC_mdata_id = 'ok'
        QC_mdata_keyword = 'Yes < 100, n= '+str(count)
    else:
        QC_mdata_id = 'good'
        QC_mdata_keyword = 'None'
        
    return QC_mdata_id, QC_mdata_keyword


def qc_Bad_Chans(infile, mad_rms, med_rms):
    """
    Checking for bad channels in the mosaic cube. 
    """

    BAD_CHAN = []

    stat_file = open(infile, 'r')
    LINES = stat_file.readlines()[2:]
    stat_file.close()

    threshold = 1.2  # value selected to be more consistent with SoFiA flagged criterion
    
#    value = med_madfm + 0.4  # Deviation from the med_madfm. Need to check with larger sample of data to decide the best value. 

    for i in range(len(LINES)):
        line = LINES[i]
        TOKS = line.split()
        chan = TOKS[0]
 #       madfm = float(TOKS[5])
        rms = float(TOKS[3])
        
#        if madfm > value:
        value = abs(rms - med_rms)
        criterion = 1.4826*threshold*mad_rms
        if value > criterion:
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


# need to see details (taken from DINGO, need to be updated for FLASH)
def qc_BeamLogs():

#    Evaluate if the synthesized beam deviate from the tolerance range
#        (+/-0.05 % of the nominal 30 arcsec which is currently set up).
#        Note that bmaj and bmin for the first few channels are always zero
#        due to barycentric correction.

#    Parameters
#    ----------
#    field: name of an interleaving field
    
    file_dir = 'SpectralCube_BeamLogs'
    basename = '/beamlog.image.restored.i.SB'+ sbid + '.cube.'+ field
    tolerance = [30 - 30 * 0.006, 30 + 30 * 0.006]

    QC_BEAMS_LABEL = []

    for i in range(0, 36):
        infile = file_dir + basename + '.beam%02d.txt' % (i)
        beamlog_file = np.loadtxt(infile)
        bmaj = beamlog_file[:,1]
        bmin = beamlog_file[:,2]
        bmaj = bmaj[bmaj > 0]
        bmin = bmin[bmin > 0]

        # check bmaj
        outliers_bmaj = (bmaj < tolerance[0]) | (bmaj > tolerance[1])

        if np.count_nonzero(outliers_bmaj) > 10:
            qc_BMAJ_label = 'fail'
        else:
            qc_BMAJ_label = 'pass'

        # check bmin
        outliers_bmin = (bmin < tolerance[0]) | (bmin > tolerance[1])
        if np.count_nonzero(outliers_bmin) > 10:
            qc_BMIN_label = 'fail'
        else:
            qc_BMIN_label = 'pass'

        # check both bmaj and bmin
        if (qc_BMAJ_label == 'pass') and (qc_BMIN_label == 'pass'):
            QC_BEAMS_LABEL.append('pass')
        else:
            QC_BEAMS_LABEL.append('fail')

    return QC_BEAMS_LABEL


def cal_beam_MADMFD(infile):
    """
    Calculating the MAD of max flux density of each beam.
    """

    data = np.loadtxt(infile)
    maxfdensity = data[:,8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 3)
    
    return mad_maxfdensity
    

def cal_beam_AvgRMS(infile):
    """
    Calculating the average RMS of each beam. 
    """
    
    data = np.loadtxt(infile)
    rms = data[:,3]
    avg_rms = round(np.mean(rms), 3)
    
    return avg_rms
    

def cal_mosaic_Stats(infile):
    """
    Calculating MAD RMS and median RMS for the mosaic cube.    
    """

    data = np.loadtxt(infile)
    rms = data[:,3]

    med_rms = np.median(rms)
    mad_rms = round(median_absolute_deviation(rms), 3)

    return mad_rms, med_rms


def cal_Theoretical_RMS(n_ant, tobs, chan_width):
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


def cal_Beam_ExpRMS(FLAGSTAT, t_rms):
    """
    Calculating the theoretical RMS of individual beam by taking into account the flagged percentage. 
    Assuming same weighting for the non-flagged data. 
    """

    BEAM_EXP_RMS = []

    for stat in FLAGSTAT:
        beam_Exp_RMS = 1/math.sqrt(float(1-stat/100)) * t_rms   # 1/sqrt(non-flagged fraction) * theoretical rms in mJy
        BEAM_EXP_RMS.append(beam_Exp_RMS)
    
    return BEAM_EXP_RMS


def cal_binnedAvg(dataArray, N):
    """
    Calculating the average value of each bin. N is the specified bin number.  
    """
    
    mean_bin = np.cumsum(dataArray, 0)[N-1::N]/float(N)
    mean_bin[1:] = mean_bin[1:] - mean_bin[:-1]
    return mean_bin



def FlagStat_plot(FLAGSTAT, n):
    """
    Plotting and visualising flagging statistics of 36 beams. 
    """

    file_dir = diagnostics_dir +'/cubestats-'+ field 
    basename = '/cubeStats-image.restored.i.SB'+ sbid +'.cube.'+ field  
    
    title = 'Flagged Fraction'
    plot_name = 'FlagStat.png'
    saved_fig = fig_dir+'/'+plot_name

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size':10}

    pylab.rcParams.update(params)

    beamXPOS, beamYPOS = BeamPosition()
    
    for i in range(36):
        bnum = n[i]
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1500, c=[FLAGSTAT[bnum]], cmap='tab20c', edgecolors='black',vmin=0, vmax=100)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

#    test with one beam - FLASH
#    plt.scatter([beamXPOS[21]], [beamYPOS[21]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 0
#    plt.scatter([beamXPOS[15]], [beamYPOS[15]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 1
#    plt.scatter([beamXPOS[14]], [beamYPOS[14]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 3

    plt.xlim(0,0.7)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    cb.set_label('Percentage')
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig)
    plt.close()

    return saved_fig, plot_name


def FlagAnt_plot(N_FLAG_ANT, n):
    """
    Plotting and visualising number of flagged (completely) antennas beam-by-beam. 
    """

    title = 'No. of 100% flagged antenna'
    plot_name = 'FlagAnt.png'
    saved_fig = fig_dir+'/'+plot_name
    
    from_list = mpl.colors.LinearSegmentedColormap.from_list

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size':10}

    pylab.rcParams.update(params)

    beamXPOS, beamYPOS = BeamPosition()

    maxflag = np.max(N_FLAG_ANT)+1
    cmap = plt.get_cmap('summer', maxflag)

    for i in range(36):
        bnum = n[i]
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1500, c=[N_FLAG_ANT[bnum]+1], cmap=cmap, edgecolors='black',vmin=0, vmax=maxflag)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

    plt.xlim(0,0.7)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    labels =np.arange(0,maxflag)
    loc = labels + .5
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig)
    plt.close()

    return saved_fig, plot_name


def Beam_ExpRMSplot(BEAM_EXPRMS, n):
    """
    Plotting and visualising expected RMS of 36 beams. 
    """

    title = 'Expected RMS'
    plot_name = 'Exp_RMS.png'
    saved_fig = fig_dir+'/'+plot_name

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size':10}

    pylab.rcParams.update(params)

    beamXPOS, beamYPOS = BeamPosition()
    VMIN = round(np.min(BEAM_EXPRMS), 3)
    VMAX = round(np.max(BEAM_EXPRMS), 3)

    for i in range(36):
        bnum = n[i]
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1500, c=[BEAM_EXPRMS[bnum]], cmap='GnBu', edgecolors='black', vmin=VMIN, vmax=VMAX)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

#    test with one beam - FLASH
#    plt.scatter([beamXPOS[21]], [beamYPOS[21]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 0
#    plt.scatter([beamXPOS[15]], [beamYPOS[15]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 1
#    plt.scatter([beamXPOS[14]], [beamYPOS[14]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 3

    plt.xlim(0,0.7)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    cb.set_label('mJy / beam')
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig)
    plt.close()

    return saved_fig, plot_name


def BeamLogs_QCplot(list_beams_id_label, n):

    params = {'axes.labelsize': 10,
              'axes.titlesize':10,
              'font.size':10}

    pylab.rcParams.update(params)
    
    legend_dict = { 'pass' : '#00FF00', 'fail' : '#CD5C5C' }
    patchList = []
    XPOS, YPOS = BeamPosition()

    plot_name = 'beamlogs_qc_SB' + sbid + '.png'
    saved_fig = fig_dir + '/' + plot_name
    
    for i in range(36):
        bnum = n[i]
        if (list_beams_id_label[bnum] =='pass'):
            color_code = '#00FF00'
        else:
            color_code = '#CD5C5C'
        
        plt.scatter([XPOS[i]], [YPOS[i]], s=1500, color=color_code, edgecolors='black')
        plt.text(XPOS[i], YPOS[i], n[i])

    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

#    test with one beam - FLASH
#    plt.scatter([XPOS[21]], [YPOS[21]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 0
#    plt.scatter([XPOS[15]], [YPOS[15]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 1
#    plt.scatter([XPOS[14]], [YPOS[14]], s=1500, c='none', edgecolors='magenta', lw=3) ##### FLASH - beam 3

    plt.legend(handles=patchList)
    plt.xlim(0,0.8)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title('Beam Log')
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def qc_NoiseRank(spread):
    """
    Evaluating the 1-percentile noise rank distribution. Sigma of the Gaussian function is used as a metric. 
    """

    sigma = 0.5 # for a Gaussian
    
    if (spread <= 2.5*sigma):
        qc_label = 'good'
    elif (2.5*sigma < spread < 3*sigma):
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
            spread = round(np.abs(coeff[2]), 3)
            ID_LABEL.append(qc_NoiseRank(spread))
#            print (infile, spread)
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

    params = {'axes.labelsize': 10,
              'axes.titlesize':10,
              'font.size':10}

    pylab.rcParams.update(params)
    
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

    params = {'axes.labelsize': 10,
              'axes.titlesize':10,
              'font.size':10}

    pylab.rcParams.update(params)

    if item == 'MADMFD':
        vmin = 0.0   # vmax = 3.0 is a conservative cut off based on M83 field. 
        vmax = 5.0
        title = 'MAD Max Flux Density'
        plot_name = 'beamStat_MADMFD.png'
        saved_fig = fig_dir+'/'+plot_name
        

    if item == 'Avg_RMS':  # this is not used 
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
    return saved_fig, plot_name



"""
def check_Cleaning ():
    
    Checking if cleaning has performed down to the defined threshold. Both major and minor cycles. 
    
    This will involve ncore specified. 
"""
    
"""
def plot(infile, x, y, c=None, yerr=None, figure=None, arrows=None, xlabel='', ylabel='')
"""


################################################################################
# Main program where it calls all the functions
################################################################################

#ignore astropy warnings 
warnings.simplefilter('ignore', AstropyWarning)   

fig_dir = 'Figures'
sbid = str(sys.argv[1]) # or, sbid=str(10250)
n = [26,25,24,23,22,21,27,10,9,8,7,20,28,11,3,1,6,19,29,12,2,0,5,18,30,13,14,15,4,17,31,32,33,34,35,16] # beam number
#n = [0] # test with one beam - FLASH

html_name = 'spectral_report_SB' + sbid + '_' + str(datetime.now().strftime("%d%b%Y")) + '.html' 

diagnostics_dir = 'diagnostics'

if not os.path.isdir(fig_dir):
    os.system('mkdir '+ fig_dir)

# Required files 

metafile = sorted(glob.glob('metadata/mslist-*txt'))[0]
metafile_science = sorted(glob.glob('metadata/mslist-scienceData*txt'))[0]
metafile_cal = sorted(glob.glob('metadata/mslist-cal*txt'))[0]
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

# Check if image cube is available. For WALLABY, this is not needed as
# the bmaj and bmin is fixed to 30"x30". 

if not os.path.isfile(fitsimage):
    bmaj = 30
    bmin = 30
else:
    bmaj, bmin = get_FitsHeader(fitsimage)

n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, total_obs_bw = get_Metadata(metafile)
askapsoft = get_Version(param)
chan_width, cfreq, nchan = get_Metadata_freq(metafile_science)
tobs_hr = round(tobs/3600.,2) # convert tobs from second to hr
chan_width_kHz = round(chan_width/1000.,3) # convert Hz to kHz
cal_sbid = metafile_cal[27:31]

cubestat_linmos_contsub = glob.glob(diagnostics_dir+ '/cubestats-' + field + '/cubeStats*linmos.contsub.txt')[0] #mosaic contsub statistic

start_freq, end_freq = get_Frequency_Range(cubestat_linmos_contsub)
freq_range = str(start_freq)+'--'+str(end_freq)
delta_freq_range = round((end_freq - start_freq), 1)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.
mad_rms, med_rms = cal_mosaic_Stats(cubestat_linmos_contsub)
n_bad_chan, mosaic_bad_chan, QC_badchan_id = qc_Bad_Chans(cubestat_linmos_contsub, mad_rms, med_rms)
QC_mdata_id, QC_mdata_keyword = qc_Missing_Data(cubestat_linmos_contsub)

first_cat = get_FIRST(ra, dec)
nvss_cat = get_NVSS(ra, dec)
#hipass_cat = get_HIPASS(ra, dec)

# Check if flagging statistic file is available. Earlier ASKAPSoft run did not include the file.
#if not os.path.isfile(flagging_file):
#    flagged_stat = 'N/A'
#else:
#    flagged_stat = get_Flagging(flagging_file)

flagging_file = sorted(glob.glob(diagnostics_dir+'/Flagging_Summaries/*averaged.ms.flagSummary')) #Flagging statistic for spectral line

FLAG_STAT, N_FLAG_ANT = [], []

if os.path.isfile(fig_dir+'/flagged_antenna.txt'):
    os.system('rm '+ fig_dir+'/flagged_antenna.txt')

for ffile in flagging_file:
    n_Rec, n_Chan, exp_count = get_Flagging_KeyValues(ffile)
    flag_stat, n_flag_ant, flag_ant_file = get_Flagging(ffile, n_Rec, n_Chan, exp_count)
    FLAG_STAT.append(flag_stat)
    N_FLAG_ANT.append(n_flag_ant)

BEAM_EXP_RMS = cal_Beam_ExpRMS(FLAG_STAT, theoretical_rms_mjy)
    
    
################################################################################  
# HTML related tasks
################################################################################

# Making thumbnail images

sizeX = 70
sizeY = 70

cube_plots = sorted(glob.glob(diagnostics_dir + '/cubestats-' + field + '/*linmos*.png')) # Mosaic statistic
beamNoise_plots = sorted(glob.glob(diagnostics_dir + '/beamNoise*.png')) # beam-by-beam statistic
beamMinMax_plots = sorted(glob.glob(diagnostics_dir +'/beamMinMax*.png')) # beam-by-beam statistic
source_spectra_plots = glob.glob('../Source_spectra/SB13293*.png') # source psectra

thumb_cubeplots = []
thumb_beamNoise = []
thumb_beamMinMax = []
thumb_source_spectra = []

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

for image in source_spectra_plots:
    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.split('/')[1]
    thumb_source_spectra.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)


beam_MADMFD_fig, MADMFD_plot = BeamStat_plot('MADMFD', n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ MADMFD_plot
make_Thumbnail(beam_MADMFD_fig, thumb_img, sizeX, sizeY, fig_dir)

# mean RMS of each beam and compares it to theoretical RMS (not taking into account flagging)
#beam_Avg_RMS_fig, AvgRMS_plot = BeamStat_plot('Avg_RMS', n)
#thumb_img = 'thumb_'+ str(sizeX) + '_'+ AvgRMS_plot
#make_Thumbnail(beam_Avg_RMS_fig, thumb_img, sizeX, sizeY, fig_dir)

# NOTE: Flagging_Summaries directory
beamExpRMS_fig, beamExpRMS_plot = Beam_ExpRMSplot(BEAM_EXP_RMS, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ beamExpRMS_plot
make_Thumbnail(beamExpRMS_fig, thumb_img, sizeX, sizeY, fig_dir)


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

list_beams_id_label = qc_BeamLogs()
BeamLogs_QCfig, BeamLogs_QCplot = BeamLogs_QCplot(list_beams_id_label, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ BeamLogs_QCplot
make_Thumbnail(BeamLogs_QCfig, thumb_img, sizeX, sizeY, fig_dir)

# NOTE: Flagging_Summaries directory
Flagged_fig, Flagged_plot = FlagStat_plot(FLAG_STAT, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ Flagged_plot
make_Thumbnail(Flagged_fig, thumb_img, sizeX, sizeY, fig_dir)

# NOTE: Flagging_Summaries directory
Flagged_ant_fig, Flagged_ant_plot = FlagAnt_plot(N_FLAG_ANT, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ Flagged_ant_plot
make_Thumbnail(Flagged_ant_fig, thumb_img, sizeX, sizeY, fig_dir)


################################################################################
# Creating html report
################################################################################

css_style = """<style>
                  body {
                        background-color: white;
                  }
                  h1 {
                      color: black;
                  }
                  h4 {
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
                    <title>ASKAP FLASH Spectral Line Data Validation Report</title>
                    """ % (css_style))

html.write("""<table width="100%" border="0">
                    <tr>
                        <th style="width:220"><img style="border:none;" src="" width="220" align="center"></th>
                        <th><h1 align="center"><font color="red">FLASH</font> Spectral Line Data Validation Report</h1>
                    	<h4 align="center">Last modified: {0} by Hyein Yoon<br>
			<font color="black">Original script for WALLABY: 24-Mar-2020 by Bi-Qing For (ICRAR/UWA)</font></h4></th>
                   	<th style="width:220"><img style="border:none;" src="" width="220" align="center"></th>
		    </tr></table>""".format(str(datetime.now().strftime("%d-%b-%Y"))))


##### Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">Notes for FLASH:
		<br>- This tool uses ASKAPsoft products. FITS-datacubes are needed for getting major and minor beam sizes only (from the header). 
		<br>- Not all data are availble, so some dummy files were used to run the script successfully.                
		<br>- 1) Combining all info from spectra + continuum  
		<br>- 2) Any other additional items to be required?
		<br><br></h4>
              </p>""")


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


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- col 1: from input by user
		<br>- col 2-8: from ./metadata/mslist-<date>*.txt
		<br>- col 9: from ./metadata/mslist-Science*.txt
		</h4>
              </p>""")
############################## 

html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Processed Image Cube</h2>
                    <table width="100%" border="1">
                    <tr>
                        <th>ASKAPsoft<br>version*</th>
                        <th>Cal SBID</th>
                        <th>Frequency Range<br>(MHz)</th>
                        <th>Central Frequency<br>(MHz)</th>
                        <th>Channel Width<br>(kHz)</th>
                        <th>Synthesised Beam<br>(arcsec x arcsec)</th>
                        <th>Beam Logs</th>
                        <th>Flagged Visibilities</th>
                        <th>Flagged Antennas</th>
                        <th>Expected RMS</th>
                    </tr>
                    <tr align="middle">
                        <td>{0}</td>
                        <td>{1}</td>
                        <td>{2}</td>
                        <td>{3}</td>
                        <td>{4}</td>
                        <td>{5}x{6}</td>
                        <td>
                        <a href="{7}" target="_blank"><img src="{8}" width="{9}" height="{10}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{11}" target="_blank"><img src="{12}" width="{13}" height="{14}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{15}" target="_blank"><img src="{16}" width="{17}" height="{18}" alt="thumbnail"></a>
                          <form action="{19}" method="get" target="_blank">
                          <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
                        </td>
                        <td>
                        <a href="{20}" target="_blank"><img src="{21}" width="{22}" height="{23}" alt="thumbnail"></a>
                        """.format(askapsoft,
                                   cal_sbid,
                                   freq_range,
                                   cfreq,
                                   chan_width_kHz,
                                   bmaj,
                                   bmin,
                                   BeamLogs_QCfig, 
                                   fig_dir+'/'+ 'thumb_' + str(sizeX) + '_' + BeamLogs_QCplot,
                                   sizeX,
                                   sizeY,
                                   Flagged_fig, 
                                   fig_dir+'/'+ 'thumb_' + str(sizeX) + '_' + Flagged_plot,
                                   sizeX,
                                   sizeY,
                                   Flagged_ant_fig, 
                                   fig_dir+'/'+ 'thumb_' + str(sizeX) + '_' + Flagged_ant_plot,
                                   sizeX,
                                   sizeY,
                                   fig_dir+'/' + str(flag_ant_file),
                                   beamExpRMS_fig, 
                                   fig_dir+'/'+ 'thumb_' + str(sizeX) + '_' + beamExpRMS_plot,
                                   sizeX,
                                   sizeY))


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- col 1: from ./slurmOutput/*.sh - if more than one version of ASKAPsoft is used for the whole reduction, the latest one is reported.
	  	<br>- col 2: from ./diagnostics/cubestats-/cubeStats*linmos.contsub.txt (mosaic contsub)
		<br>- col 3-4: from ./metadata/mslist-Science*.txt
		<br>- col 5: from FITS-datacube (CURRENT VERSION: continuum subtracted beam00 cube - Nov 22 ver.; too large beam size? depending on robust parameter?)
		<br>- col 6: from ./SpectralCube_BeamLogs/beamlogs*.txt
		<br>- col 6: Bi-qing's notes: Evaluating each channel of each beam if ASKAPSoft fails to synthesize the beam, bmaj and bmin to 30 arcsec. bmaj and bmin for the first few channels are always zero.
		<br>- col 7: from ./flagSummary/*.flagSummary
		<br>- col 8: from ./flagSummary/*.flagSummary (flagged fraction) + theoretical rms estimation (based on input values)		
		</h4> 
              </p>""")


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
                               beamNoise_plots[1],
                               fig_dir+'/'+ thumb_beamNoise[1],
                               sizeX,
                               sizeY,
                               beamNoise_plots[0],
                               fig_dir+'/'+ thumb_beamNoise[0],
                               sizeX,
                               sizeY,
                               beamNoise_plots[2],
                               fig_dir+'/'+ thumb_beamNoise[2],
                               sizeX,
                               sizeY))


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- col 1-3: from beamMinMax Plots 
		<br>- why one percentile?</h4>
              </p>""")


html.write("""</td>
                    </tr>
                    </table>
                    <table width="100%" border="1">
                    <tr>
                    <th colspan="4">Continuum Subtracted Beam Cube</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>MAD Max Flux Density</p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    <br><p>1-percentile noise rank</p>
                    """.format(fig_dir+'/'+'beamStat_MADMFD.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_MADMFD.png',
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


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- col 1: from beamMinMax Plots 
	 	<br>- col 2: from CubeStat*contsub.txt</h4>
              </p>""")


html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Mosaic Statistics</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Image Cube</th>          
                    <th>Continuum Subtracted Cube</th>
                    <th>Residual Cube</th>
                    <th>Number of Bad Channel</th>
                    <th>Missing Data <br>(Channel)</th>
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
                    <td id='{12}'>{13}
                    <form action="{14}" method="get" target="_blank">
                     <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
                    </form>
                    <td id='{15}'>{16}
                    """.format(cube_plots[1],
                               fig_dir+'/' + thumb_cubeplots[1],
                               sizeX,
                               sizeY,
                               cube_plots[0],
                               fig_dir+'/'+ thumb_cubeplots[0],
                               sizeX,
                               sizeY,
                               cube_plots[2],
                               fig_dir+'/'+ thumb_cubeplots[2],
                               sizeX,
                               sizeY,
                               QC_badchan_id,
                               n_bad_chan,
                               fig_dir+'/' + mosaic_bad_chan,
                               QC_mdata_id,
                               QC_mdata_keyword))


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- col 1-3: from cubePlots
	      <br>- col 4: from CubeStat*contsub.txt</h4>
              </p>""")


html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Source and Noise Spectra from five bright components</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Component 01a</th>          
                    <th>Component 01b</th>
                    <th>Component 02a</th>
                    <th>Component 03a</th>
                    <th>Component 03b</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="../Source_spectra/SB13293_component_1a.png" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p></p>
		    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p></p>
		    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p></p>
		    </td>
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p></p>
                    </td>
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    <br><p></p>
	            </td>
<tr align="middle">
                    <td>
                    <a href="../Source_spectra/SB13293_component_1a.png" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>10/33 chunks > 5-sigma</p>
		    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>13/33 chunks > 5-sigma</p>
		    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>13/33 chunks > 5-sigma</p>
		    </td>
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p>9/33 chunks > 5-sigma</p>
                    </td>
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    <br><p>11/33 chunks > 5-sigma</p>
	            </td>

                    """.format("SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_1a.png",
                               sizeX,
                               sizeY,
                               "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_1b.png",
                               sizeX,
                               sizeY,
                               "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_2a.png",
                               sizeX,
                               sizeY,
			        "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_2a.png",
                               sizeX,
                               sizeY,
			       "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_2a.png",
                               sizeX,
                               sizeY,
			       "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_2a.png",
                               sizeX,
                               sizeY,
			       "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_2a.png",
                               sizeX,
                               sizeY,
			       "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_2a.png",
                               sizeX,
                               sizeY,

                               "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_3a.png",
                               sizeX,
                               sizeY,
                               "SB13293_component_1a.png",
                               './Source_spectra/' + "SB13293_component_3b.png",
                               sizeX,
                               sizeY))                    



############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- Spectra toward five brightest components
              <br>- Deviation from noise spectra (9 MHz chunks)
	      </h4>
              </p>""")


# noise spectra offset
html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Median noise flux density - noise Spectra</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Low frequency (first 5,000 channels)</th>          
                    <th>High frequency (last 5,000 channels)</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>199 component (outside 3.2 deg)</p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>139 component (outside 3.2 deg)</p>
                    </td>
                    <tr align="middle">
                    <td>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    <br><p>RA offset (red points: outside 3.2 deg)</p>
                    </td>
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p>RA offset (red points: outside 3.2 deg)</p>
                    </td>
                    <tr align="middle">
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    <br><p>DEC offset (red points: outside 3.2 deg)</p>
                    </td>
                    <td>
                    <a href="{20}" target="_blank"><img src="{21}" width="{22}" height="{23}" alt="thumbnail"></a>
                    <br><p>DEC offset (red points: outside 3.2 deg)</p>
                    """.format("noisecolor_vs_radecoffset_sb13375_lowfreq.png",
                               fig_dir+'/'+ "noisecolor_vs_radecoffset_sb13375_lowfreq.png",
                               sizeX,
                               sizeY,
                               "noisecolor_vs_radecoffset_sb13375_highfreq.png",
                               fig_dir+'/'+ "noisecolor_vs_radecoffset_sb13375_highfreq.png",
                               sizeX,
                               sizeY,
                               "noise_vs_raoffset_outliers_sb13375_lowfreq.png",
                               fig_dir+'/'+ "noise_vs_raoffset_outliers_sb13375_lowfreq.png",
                               sizeX,
                               sizeY,
                               "noise_vs_raoffset_outliers_sb13375_highfreq.png",
                               fig_dir+'/'+ "noise_vs_raoffset_outliers_sb13375_highfreq.png",
                               sizeX,
                               sizeY,
                               "noise_vs_decoffset_outliers_sb13375_lowfreq.png",
                               fig_dir+'/'+ "noise_vs_decoffset_outliers_sb13375_lowfreq.png",
                               sizeX,
                               sizeY,
                               "noise_vs_decoffset_outliers_sb13375_highfreq.png",
                               fig_dir+'/'+ "noise_vs_decoffset_outliers_sb13375_highfreq.png",
                               sizeX,
                               sizeY))

############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- Mean noise flux density - noise spectra 
		<br>- stable out to 3.2 degree</h4>
              </p>""")


html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Continuum - comparison with NVSS</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Continuum image</th>          
                    <th>Statistics</th>
                    <th>RA/DEC offset</th>
                    <th>Flux comparison</th>
                    <th>Flux vs distance from image centre</th>
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
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    """.format("nvss_askap_components_sb10250.png",
                               fig_dir+'/' + "nvss_askap_components_sb10250.png",
                               sizeX,
                               sizeY,
                               "continuum_statitiscs.png",
                               fig_dir+'/'+ "continuum_statitiscs.png",
                               sizeX,
                               sizeY,
                               "crossmatch_ra_dec_offset_sb10849.png",
                               fig_dir+'/'+ "crossmatch_ra_dec_offset_sb10849.png",
                               sizeX,
                               sizeY,
                               "crossmatch_flux_ratio_sb10849.png",
                               fig_dir+'/'+ "crossmatch_flux_ratio_sb10849.png",
                               sizeX,
                               sizeY,
                               "crossmatch_flux_vs_image_center_sb10849.png",
                               fig_dir+'/' + "crossmatch_flux_vs_image_center_sb10849.png",
 			       sizeX,
                               sizeY))


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- col 1: continuum image + selavy bright componenets
	      <br>- col 2: size & flux histogram
	      <br>- col 3: RA/DEC offset (comparison with NVSS)
	      <br>- col 4: flux difference (comparison with NVSS)
              <br>- col 5: primary beam correction check (comparison with NVSS)</h4>
              </p>""")



##############################  for a quick test - remove below
html.write("""</td> 
              </tr>
              </table>
              <h2 align="middle">FIRST sources within 6 x 6 sq degree</h2>
              <table>
              <tr align="middle">
              <td id="cell">
              <form action="{0}" method="get" target="_blank">
                 <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
              </form>
              """.format(fig_dir+'/' + first_cat))
############################## 


##############################  for a quick test - remove below
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- data from Vizier FIRST (2014Dec17; Helfand+ 2015)
              <br>- a resolution of 5 arcsec</h4>
              </p>""")

html.write("""</td> 
              </tr>
              </table>
              <h2 align="middle">NVSS sources within 6 x 6 sq degree</h2>
              <table>
              <tr align="middle">
              <td id="cell">
              <form action="{0}" method="get" target="_blank">
                 <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
              </form>
              """.format(fig_dir+'/' + nvss_cat))
############################## 


############################## Notes for FLASH
html.write("""
              </td>
              </tr>
              </table>
                                             
              <p align=left>
              <h4 align="left">- data from Vizier NVSS (Condon+ 1998)
              <br>- a resolution of 45 arcsec</h4>
              </p>""")



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
              <a href="mailto:hyein.yoon@sydney.edu.au">Hyein Yoon</a></i>
              </p>

              </body>
              </html>""".format(str(datetime.now())))

html.close()
print ("Spectral line validation report written to '{0}'.".format(html_name))
