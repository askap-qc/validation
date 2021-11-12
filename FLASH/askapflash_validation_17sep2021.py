################################################################################
# ASKAP-FLASH data validation script
################################################################################

################################################################################
#
# This script generates an ASKAP-FLASH spectral line cube HTML report. 
# Compatibility: Python version 3.x
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- metadata/mslist-cal*txt
# -- metadata/mslist-*101.txt
# -- slurmOutput/*sh
# -- diagnostics/cubestats-<field>/*txt
# -- diagnostics/*png
# -- diagnostics/Flagging_Summaries/*SL.ms.flagSummary 
# -- SpectralCube_BeamLogs/*.txt
# -- image.restored.i.SB<SBID>.cube.contsub.fits (to read bmaj, bmin from header)
#
# Output files: all saved in Figures directory. 
#
# To run:
#
# python <script name> -s <SBID> -c <cal_SBID>
# e.g. python askapflash_validation.py -s 13293 -c 13293
#
# The original script was written for ASKAP-WALLABY data validation by Bi-Qing For 
# Version: 24 March 2020
#
# Updated for ASKAP-FLASH data validation by Hyein Yoon
# Email: hyein.yoon [at] sydney.edu.au
# Last update: 17 September 2021
#
################################################################################


################################################################################
# Modules for the script
################################################################################

import os
import os.path
import sys
import glob
import math
import subprocess 
import numpy as np
import PIL 
import warnings
import astropy.coordinates as coord
import astropy.units as u

from numpy import inf
from datetime import datetime
from PIL import Image

from astropy.io import fits
from astropy.io.fits import getheader
from astropy.stats import median_absolute_deviation
from astropy.utils.exceptions import AstropyWarning

from scipy.optimize import OptimizeWarning
from scipy.constants import k as k_B
from scipy.stats import iqr
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

from astroquery.vizier import Vizier
from argparse import ArgumentParser

import matplotlib as mpl
mpl.use('Agg')  # to avoid matplotlib using the Xwindows backend

import matplotlib.pyplot as plt 
import matplotlib.pylab as pylab
import matplotlib.patches as mpatches


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

    return A*np.exp(-(x-mu)**2.0/(2.0*sigma**2.0))


def get_FitsHeader(fitsimage):
    """
    Getting basic information of processed mosaic.
    """

    hdu = getheader(fitsimage)
    bmaj = round(float(hdu['BMAJ'])*3600., 2)  # arcsec
    bmin = round(float(hdu['BMIN'])*3600., 2)  # arcsec

    return bmaj, bmin


def get_NVSS(ra, dec):
    """
    Getting the NVSS sources within the 6x6 sq deg field through VizieR.
    NVSS: Condon+ 1998.
    """

    from astroquery.vizier import Vizier # for installation: conda install -c astropy astroquery
 
    print ("Retrieving NVSS sources from Vizier. Depending on server connection, this might take a while....")

    catalogue='VIII/65/nvss'

    v = Vizier(columns=['NVSS', '_RAJ2000', '_DEJ2000', 'S1.4', 'MajAxis', 'MinAxis'], catalog = catalogue, timeout=10000)
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
    if nvss_result.keys()==[catalogue]:
#      print (nvss_result[catalogue], file=open(fig_dir + '/' + nvss_cat,'w'))
       nvss_result['VIII/65/nvss'].write(fig_dir + '/' + nvss_cat, format='ascii.fixed_width', delimiter=' ')
    else:
        print ("No sources found", file=open(fig_dir + '/' + nvss_cat,'w'))

    return nvss_cat


def get_Version(param):
    """
    Getting the latest ASKAPsoft version that is used for the data reduction.
    """

    line = subprocess.check_output(['grep', 'Processed with ASKAPsoft', param])
    str_line = line.decode('utf-8')

    newline = str_line.splitlines()[-1] #Get the most recent, in case there is more than one instance
    askapsoft = newline.split()[-1]

    return askapsoft
        

def get_Flagging_KeyValues(flagging_file):
    """
    Getting Flagging Key Values. 
    """

    flag_infile = open(flagging_file, 'r')
    LINES = flag_infile.readlines()[:10]
    flag_infile.close()
    
    N_Rec = 'nRec'  # Total number of spectra feeds into the synthesis image. This is not always constant so grab the value beam-by-beam.
    N_Chan = 'nChan'  # Total number of channel

    # Searching for keywords in the file   
    for i in range(len(LINES)):
        line = LINES[i]
        if line.find(N_Rec) >=0:
            TOKS = line.split()
            n_Rec = float(TOKS[2])
        if line.find(N_Chan) >=0:
            TOKS = line.split()
            n_Chan = float(TOKS[2])

    exp_count = n_Rec*35  #counting antenna number from zero based on the recorded data
    
    return n_Rec, n_Chan, exp_count


def get_Flagging(flagging_file, n_Rec, nChan, exp_count):
    """
    Getting flagging statistics and finding out beam-by-beam antenna based (completely) flagging. 
    """

    line = subprocess.check_output(['grep','Flagged', flagging_file]) #grab the summary line
    str_line = line.decode('utf-8')
    TOKS = str_line.split()
    total_flagged_pct = float(TOKS[-2]) #data + autocorrelation
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
                    if len(line.split())>=7:  # avoid new channel-wise summaries at end of flagSummary file
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
        ffile.write(flagging_file[-30:-24])
        ffile.write('\n')
        for item in ANT_NAME:
            ffile.write(item)
            ffile.write('\n')
    else:
        ffile.write(flagging_file[-30:-24])
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
    
    # FLASH: bandpass smoothing affected 18 sbids - processed prior to Mar 2021
    bs_issue_sbid=[13285,13283,13334,13335,13336,13273,13297,13296,13278,15873,13271,10849,13293,11068,13306,13294,15208,15209]

    bs_issue_id="No"
    for j in range(len(bs_issue_sbid)):
        if bs_issue_sbid[j] == float(sbid):
            bs_issue_id="Yes"
            break

    return n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, total_obs_bw, bs_issue_id


def get_Metadata_freq(metafile_science):
    """
    Getting basic information on observed field (one field only). 
    Frequency and channel related items. 
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
    Frequency range and total channel of the mosaic.
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
    Making thumbnail images.
    """
    
    size = sizeX, sizeY
    im = Image.open(image)
    im.thumbnail(size)
    im.save(dir+ '/' + thumb_img)

    return thumb_img


def qc_Missing_Data(infile):
    """
    Checking for missing data channels in the mosaic cube.
    Printing the list of missing data channels. 
    """

    data = open(infile, 'r')
    LINES = data.readlines()[2:]

    MDATA_CHAN = []

    for i in range(len(LINES)):
        line = LINES[i]
        TOKS = line.split()
        chan = TOKS[0]
        maxfdensity = float(TOKS[8])

        if maxfdensity==0:
            MDATA_CHAN.append(chan)

    if MDATA_CHAN == []:
        MDATA_CHAN.append('none')
        QC_mdata_chan_id = 'good'
        QC_mdata_chan_keyword = 'n= '+str(len(MDATA_CHAN))
    else:
        QC_mdata_chan_id = 'bad'
        QC_mdata_chan_keyword = 'n= '+str(len(MDATA_CHAN))

    mosaic_mdata_chan = 'mosaic_missing_data.txt'
    print ('\n'.join(MDATA_CHAN), file=open(fig_dir + '/' + mosaic_mdata_chan,'w'))

    n_mdata_chan = len(MDATA_CHAN)

    # Check if number of bad channel recorded is 1. If yes, check if is it a none keyword.
    # If yes, number of bad channel should be 0.
    
    if n_mdata_chan == 1:
        with open(fig_dir + '/' + mosaic_mdata_chan) as f:
            if 'none' in f.read():
                n_mdata_chan = 0
                print ('yes')
        
    return  n_mdata_chan, mosaic_mdata_chan, QC_mdata_chan_id, QC_mdata_chan_keyword


def qc_Bad_Chans(infile, mad_rms, med_rms):
    """
    Checking for bad data channels in the mosaic cube. 
    Printing the list of bad data channels. 
    """

    data = np.loadtxt(infile)
    data_rms = data[:,3]
    data_median_rms=np.median(data_rms)
    data_std_rms=np.std(data_rms)
    threshold = np.std(data_rms)*3.0  # threshold = 3 * standard deviation from med_rms

    BAD_CHAN = []

    stat_file = open(infile, 'r')
    LINES = stat_file.readlines()[2:]
    stat_file.close()

    for i in range(len(LINES)):
        line = LINES[i]
        TOKS = line.split()
        chan = TOKS[0]
        freq = TOKS[1]
        rms = float(TOKS[3])
        madfm = float(TOKS[5])
        
        if abs(rms-data_median_rms) > threshold:
            BAD_CHAN.append(chan)
        elif rms == 0:
            BAD_CHAN.append(chan)

    if BAD_CHAN == []:
        BAD_CHAN.append('none')
        QC_badchan_id = 'good'
        QC_badchan_keyword = 'n= '+str(len(BAD_CHAN))
    else:
        QC_badchan_id = 'bad'
        QC_badchan_keyword = 'n= '+str(len(BAD_CHAN))


    mosaic_bad_chan = 'mosaic_bad_chan.txt'
    print ('\n'.join(BAD_CHAN), file=open(fig_dir + '/' + mosaic_bad_chan,'w'))

    n_bad_chan = len(BAD_CHAN)

    # Check if number of bad channel recorded is 1. If yes, check if is it a none keyword.
    # If yes, number of bad channel should be 0.
    
    if n_bad_chan == 1:
        with open(fig_dir + '/' + mosaic_bad_chan) as f:
            if 'none' in f.read():
                n_bad_chan = 0
                print ('yes')
    
    return n_bad_chan, mosaic_bad_chan, QC_badchan_id, QC_badchan_keyword


def qc_BeamLogs():
    """
     Threshold - 30% deviation from median bmaj, bmin of each beam. 
     For FLASH, beam sizes change greatly across the full 288 MHz band.
     Checking with some plots - e.g. bmaj or bmin vs chan is recommended.
     Note that bmaj and bmin for the first few channels are always zero due to barycentric correction.
    """
  
    file_dir = 'SpectralCube_BeamLogs'
    basename = '/beamlog.image.restored.' + imagebase + '.' + field

    # use different basename - latest sbids
    if not glob.glob(file_dir + basename +'*.txt'):  
        basename = '/beamlog.image.restored.' + 'i.' + field + '.SB' + sbid + '.cube.' + field  

    QC_BEAMS_LABEL = []
    MEDIAN_BMAJ = []
    MEDIAN_BMIN = []

    for i in range(0, 36):
        infile = file_dir + basename + '.beam%02d.txt' % (i)
        if os.path.isfile(infile):
            beamlog_file = np.loadtxt(infile)
            bmaj = beamlog_file[:,1]
            bmin = beamlog_file[:,2]
            bmaj = bmaj[bmaj > 0]
            bmin = bmin[bmin > 0]
            bmaj_median=np.median(bmaj)
            bmin_median=np.median(bmin)
            tolerance_maj=[bmaj_median - bmaj_median*0.3, bmaj_median + bmaj_median*0.3]
            tolerance_min=[bmin_median - bmin_median*0.3, bmin_median + bmin_median*0.3]
            MEDIAN_BMAJ.append(bmaj_median)
            MEDIAN_BMIN.append(bmin_median)

        # check bmaj
            outliers_bmaj = (bmaj < tolerance_maj[0]) | (bmaj > tolerance_maj[1])

            if np.count_nonzero(outliers_bmaj) > 50:
                qc_BMAJ_label = 'fail'
            else:
                qc_BMAJ_label = 'pass'

        # check bmin
            outliers_bmin = (bmin < tolerance_min[0]) | (bmin > tolerance_min[1])

            if np.count_nonzero(outliers_bmin) > 50:
                qc_BMIN_label = 'fail'
            else:
                qc_BMIN_label = 'pass'

        # check both bmaj and bmin
            if (qc_BMAJ_label == 'pass') and (qc_BMIN_label == 'pass'):
                QC_BEAMS_LABEL.append('pass')
            else:
                QC_BEAMS_LABEL.append('fail')

        # no beamlogs text files - missing beams
        else:
                QC_BEAMS_LABEL.append('missing')
                MEDIAN_BMAJ.append(0)
                MEDIAN_BMIN.append(0)

    return QC_BEAMS_LABEL, MEDIAN_BMAJ, MEDIAN_BMIN


def cal_beam_MADMFD(infile):
    """
    Calculating the MAD of max flux density of each beam.
    """

    data = np.loadtxt(infile)
    maxfdensity = data[:,8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 3)
    
    return mad_maxfdensity
    

def cal_beam_MeanRMS(infile):
    """
    Calculating the average RMS of each beam. 
    """
    
    data = np.loadtxt(infile)
    rms = data[:,3]
    mean_rms = round(np.mean(rms), 3)
    
    return mean_rms


def cal_beam_MedianRMS(infile):
    """
    Calculating the median RMS of each beam. 
    """
    
    data = np.loadtxt(infile)
    rms = data[:,3]
    median_rms = round(np.median(rms), 3)
    
    return median_rms


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
    Calculating the theoretical rms noise for ASKAP. 
    Assuming natural weighting and not taking into account fraction of flagged data. 
    """

    tsys = 60       # K # tsys/eta ~ 85K (Band 1 average, Hotan+ 2021 Fig.22)
    antdiam = 12    # metre
    aper_eff = 0.7  # aperture efficiency # 0.7 - 0.75
    coreff = 0.8    # correlator efficiency
    npol = 2.0      # Number of polarisation, npol = 2 for images in Stokes I, Q, U, or V
    
    anteff = np.pi*(antdiam/2.0)**2.0*aper_eff
    SEFD = 2.0*k_B*1e26*tsys/anteff   
    rms_jy = SEFD/(coreff*np.sqrt(npol*n_ant*(n_ant-1.0)*chan_width*tobs))

    return rms_jy


def cal_Beam_ExpRMS(FLAGSTAT, t_rms):
    """
    Calculating the theoretical RMS of individual beam by taking into account the flagged percentage. 
    Assuming same weighting for the non-flagged data. 
    """

    BEAM_EXP_RMS = []

    for stat in FLAGSTAT:
        # 1/sqrt(non-flagged fraction) * theoretical rms in mJy
        beam_Exp_RMS = 1/math.sqrt(float(1-stat/100)) * t_rms   
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
    basename = '/cubeStats-image.restored.' + imagebase + '.' + field  

    # use different basename - latest sbids
    if not glob.glob(file_dir + basename +'*.txt'):
        basename = '/cubeStats-image.restored.' + 'i.' + field + '.SB' + sbid + '.cube.' + field  
    
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
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1400, c=[FLAGSTAT[bnum]], cmap='tab20c', edgecolors='black',vmin=0, vmax=100)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

    plt.xlim(0,0.7)
    plt.ylim(0,1.4) 

    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    cb.set_label('Percentage')
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig, bbox_inches='tight') 
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
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1400, c=[N_FLAG_ANT[bnum]], cmap=cmap, edgecolors='black',vmin=0, vmax=maxflag)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

    plt.xlim(0,0.7)
    plt.ylim(0,1.4) 

    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    labels =np.arange(0,maxflag)
    loc = labels + .5
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig, bbox_inches='tight') 
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
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1400, c=[BEAM_EXPRMS[bnum]], cmap='GnBu', edgecolors='black', vmin=VMIN, vmax=VMAX)
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

    plt.xlim(0,0.7)
    plt.ylim(0,1.4) 
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    cb.set_label('mJy / beam')
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig, bbox_inches='tight') 
    plt.close()

    return saved_fig, plot_name


def BeamLogs_QCplot(list_beams_id_label, n):
    """
    Plotting and visualising beamlogs of 36 beams. 
    """

    params = {'axes.labelsize': 10,
              'axes.titlesize':10,
              'font.size':10}

    pylab.rcParams.update(params)
    
    legend_dict = { 'pass' : '#00FF00', 'fail_partial' : '#CD5C5C', 'missing' : '#E4E4E4'}
    patchList = []
    XPOS, YPOS = BeamPosition()

    plot_name = 'beamlogs_qc_SB' + sbid + '.png'
    saved_fig = fig_dir + '/' + plot_name
    
    for i in range(36):
        bnum = n[i]
        if (list_beams_id_label[bnum] =='pass'):  #index 0 corresponds to beam 26 
            color_code = '#00FF00'
        elif (list_beams_id_label[bnum] =='missing'):  # for missing beams
            color_code = '#E4E4E4'
        else:
            color_code = '#CD5C5C'

        plt.scatter([XPOS[i]], [YPOS[i]], s=1400, color=color_code, edgecolors='black')
        plt.text(XPOS[i], YPOS[i], n[i])

    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    plt.legend(handles=patchList,loc='best')
    plt.xlim(0.01,0.91) 
    plt.ylim(0,1.39) 

    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title('Beam Log')
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def qc_NoiseRank(spread):
    """
    Evaluating the 1-percentile noise rank distribution. 
    Sigma of the Gaussian function is used as a metric. 
    """

    sigma = 0.5 # Gaussian distribution (default value, sigma=0.5)
    
    if (spread <= 3*sigma):
        qc_label = 'good'
    elif (3*sigma < spread < 6*sigma):
        qc_label = 'ok'
    else:
        qc_label = 'uncertain'

    return qc_label


def NoiseRank_histplot(nchan):
    """
    Plotting and visualising 1-percentile noise level histogram of 36 beams. 
    """
    
    ID_LABEL = []
    plot_name = 'beam_1pctile_hist_SB'+ sbid + '.png'
    saved_fig = fig_dir + '/' + plot_name
    file_dir = diagnostics_dir +'/cubestats-'+ field 
    basename = '/cubeStats-image.restored.' + imagebase + '.' + field

    # use different basename - latest sbids
    if not glob.glob(file_dir + basename +'*.txt'):
        basename = '/cubeStats-image.restored.' + 'i.' + field + '.SB' + sbid + '.cube.' + field  

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

            axs[i].bar(xcenter, N)
            axs[i].plot(xcenter,gauss(xcenter,*coeff),'r-',lw=1)    
            axs[i].set_xlim(xmin_val-3, xmax_val+3)
            axs[i].set_ylim(0, ymax_val+3)
            axs[i].title.set_text('Beam%02d' %(i))

    plt.tight_layout()
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name, ID_LABEL


def NoiseRank_QCplot(list_id_label, n):
    """
    Plotting and visualising 1-percentile noise rank of 36 beams. 
    """

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
        
        plt.scatter([XPOS[i]], [YPOS[i]], s=1400, color=color_code, edgecolors='black')
        plt.text(XPOS[i], YPOS[i], n[i])

    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    plt.legend(handles=patchList)
    plt.xlim(-0.04,0.91) 
    plt.ylim(0,1.39) 
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title('1-percentile noise rank')
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def BeamStat_plot(item, n):
    """
    Plotting and visualising statistics of 36 beams. 
    """

    file_dir = diagnostics_dir +'/cubestats-'+ field 
    basename = '/cubeStats-image.restored.' + imagebase + '.' + field  

    # use different basename - latest sbids
    if not glob.glob(file_dir + basename +'*.txt'):
        basename = '/cubeStats-image.restored.' + 'i.' + field + '.SB' + sbid + '.cube.' + field  

    params = {'axes.labelsize': 10,
              'axes.titlesize':10,
              'font.size':10}

    pylab.rcParams.update(params)

    if item == 'Max': 
        title = 'MAD Max Flux Density'
        plot_name = 'beamStat_MADMFD.png'
        saved_fig = fig_dir+'/'+plot_name
        
    if item == 'Mean':           
        title = 'Mean RMS'
        plot_name = 'beamStat_MeanRMS.png'
        saved_fig = fig_dir+'/'+plot_name

    if item == 'Median':   
        title = 'Median RMS'
        plot_name = 'beamStat_MedianRMS.png'
        saved_fig = fig_dir+'/'+plot_name
    
    beamXPOS, beamYPOS = BeamPosition()
    
    beamstat_all = [] 
    beamstat_only_data = []

    for i in range(36):
        bnum = n[i]
        infile = file_dir + basename +'.beam%02d.contsub.txt'%(bnum)

        if os.path.isfile(infile):
            if item == 'Max':
                beamstat = cal_beam_MADMFD(infile)
            elif item == 'Mean':	               
                beamstat = cal_beam_MeanRMS(infile)
            elif item == 'Median':
                beamstat = cal_beam_MedianRMS(infile)
            beamstat_only_data.append(beamstat)  # for vmin and vmax (excluding missing beam)

        else:  # missing beams
            print('Failed to read cubeStats field for beams')
            beamstat = 99.0

        beamstat_all.append(beamstat)

        # vmin=min(beamstat_only_data)
        # vmax=max(beamstat_only_data)
 
        plt.scatter([beamXPOS[i]], [beamYPOS[i]], s=1400, c=[beamstat_all[i]], cmap='RdYlGn_r', edgecolors='black', vmin=0, vmax=max(beamstat_only_data))
        plt.text(beamXPOS[i], beamYPOS[i], n[i])

    plt.xlim(0,0.7)
    plt.ylim(0,1.4)
    plt.tick_params(axis='both',which='both', bottom=False,top=False,right=False,left=False,labelbottom=False, labelleft=False)
    plt.title(title)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=10)
    cb.set_label('mJy / beam')
    plt.savefig(saved_fig,bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


################################################################################
# Main program where it calls all the functions
################################################################################


# ignore astropy warnings 
warnings.simplefilter('ignore', AstropyWarning)   

parser = ArgumentParser(description='Run FLASH validation and produce an HTML report')
parser.add_argument('-s','--sbid', dest='sbid',required='true',help='Science SBID',type=str)
parser.add_argument('-c','--cal_sbid', dest='cal_sbid',required='true',help='Calibrator SBID',type=str)
parser.add_argument('-i','--imagebase', dest='imagebase',default='i.SB%s.cube',help='Base string for images [default=%default]',type=str)
options = parser.parse_args()

fig_dir = 'Figures'
sbid = options.sbid
cal_sbid = options.cal_sbid

# beam number
n = [26,25,24,23,22,21,27,10,9,8,7,20,28,11,3,1,6,19,29,12,2,0,5,18,30,13,14,15,4,17,31,32,33,34,35,16] 

html_name = 'index.html' 

diagnostics_dir = 'diagnostics'

imagebase = options.imagebase
imagebase = imagebase.replace('%s',sbid)

if not os.path.isdir(fig_dir):
    os.system('mkdir '+ fig_dir)

# Required files 

metafile = sorted(glob.glob('metadata/mslist-*txt'))[0]
metafile_science = sorted(glob.glob('metadata/mslist-scienceData*txt'))[0]
param_file = sorted(glob.glob('slurmOutput/*.sh'))
fitsimage = ('image.restored.' + imagebase + '.contsub.fits')

# Check if there is more than one parameter input .sh file in the slurmOutput directory.
# If it does, select the latest one.
# If more than one version is used. Reporting the latest version of ASKAPSoft for final data reduction. 

if len(param_file) >=1:
    index = len(param_file)
    param = param_file[index-1]
else:
    param = param_file[0]


n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, total_obs_bw, bs_issue_id = get_Metadata(metafile)
askapsoft = get_Version(param)
chan_width, cfreq, nchan = get_Metadata_freq(metafile_science)
tobs_hr = round(tobs/3600.,2) # convert tobs from second to hr
chan_width_kHz = round(chan_width/1000.,3) # convert Hz to kHz

# Beamsize - check if image cube is available. 
# if not, calculate median of bmaj and bmin from beamlogs

list_beams_id_label, median_bmaj, median_bmin = qc_BeamLogs()

if not os.path.isfile(fitsimage):
    header_bmaj = round(np.median(median_bmaj),2)
    header_bmin = round(np.median(median_bmin),2)
else:
    header_bmaj, header_bmin = get_FitsHeader(fitsimage)

cubestat_linmos_contsub = glob.glob(diagnostics_dir+ '/cubestats-' + field + '/cubeStats*linmos.contsub.txt')[0] #mosaic contsub statistic

start_freq, end_freq = get_Frequency_Range(cubestat_linmos_contsub)
freq_range = str(start_freq)+'--'+str(end_freq)
delta_freq_range = round((end_freq - start_freq), 1)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.
mad_rms, med_rms = cal_mosaic_Stats(cubestat_linmos_contsub)
n_bad_chan, mosaic_bad_chan, QC_badchan_id, QC_badchan_keyword = qc_Bad_Chans(cubestat_linmos_contsub, mad_rms, med_rms)
n_mdata_chan, mosaic_mdata_chan, QC_mdata_chan_id, QC_mdata_chan_keyword = qc_Missing_Data(cubestat_linmos_contsub)

# first_cat = get_FIRST(ra, dec)
nvss_cat = get_NVSS(ra, dec)

#Flagging statistic for spectral line
flagging_file = sorted(glob.glob(diagnostics_dir+'/Flagging_Summaries/*averaged.ms.flagSummary')) 

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

cube_plots = sorted(glob.glob(diagnostics_dir + '/cubestats-' + field + '/*.cube*linmos*.png'))  #Mosaic statistic
beamNoise_plots = sorted(glob.glob(diagnostics_dir + '/beamNoise*.cube*.png')) #beam-by-beam statistic
beamMinMax_plots = sorted(glob.glob(diagnostics_dir +'/beamMinMax*.cube*.png')) #beam-by-beam statistic
#source_spectra_plots = glob.glob('../Source_spectra/SB13293*.png') # source spectra

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

#for image in source_spectra_plots:
#    thumb_img = 'thumb_'+ str(sizeX) + '_'+ image.split('/')[1]
#    thumb_source_spectra.append(thumb_img)
#    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)


beam_MADMFD_fig, MADMFD_plot = BeamStat_plot('Max', n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ MADMFD_plot
make_Thumbnail(beam_MADMFD_fig, thumb_img, sizeX, sizeY, fig_dir)

# mean RMS of each beam and compares it to theoretical RMS (not taking into account flagging)
beam_Mean_RMS_fig, MeanRMS_plot = BeamStat_plot('Mean', n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ MeanRMS_plot
make_Thumbnail(beam_Mean_RMS_fig, thumb_img, sizeX, sizeY, fig_dir)

# median RMS of each beam and compares it to theoretical RMS (not taking into account flagging)
beam_Median_RMS_fig, MedianRMS_plot = BeamStat_plot('Median', n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ MedianRMS_plot
make_Thumbnail(beam_Median_RMS_fig, thumb_img, sizeX, sizeY, fig_dir)

# NOTE: careful - the location of Flagging_Summaries directory
beamExpRMS_fig, beamExpRMS_plot = Beam_ExpRMSplot(BEAM_EXP_RMS, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ beamExpRMS_plot
make_Thumbnail(beamExpRMS_fig, thumb_img, sizeX, sizeY, fig_dir)

beam_1pctile_histfig, onepctile_plot, list_id_label = NoiseRank_histplot(float(nchan))
thumb_img = 'thumb_' + str(sizeX) + '_' + onepctile_plot
make_Thumbnail(beam_1pctile_histfig, thumb_img, sizeX, sizeY, fig_dir)

beam_1pctile_QCfig, onepctile_QCplot = NoiseRank_QCplot(list_id_label, n)
thumb_img = 'thumb_' + str(sizeX) + '_' + onepctile_QCplot
make_Thumbnail(beam_1pctile_QCfig, thumb_img, sizeX, sizeY, fig_dir)

list_beams_id_label, median_bmaj, median_bmin = qc_BeamLogs()
BeamLogs_QCfig, BeamLogs_QCplot = BeamLogs_QCplot(list_beams_id_label, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ BeamLogs_QCplot
make_Thumbnail(BeamLogs_QCfig, thumb_img, sizeX, sizeY, fig_dir)

# NOTE: careful - the location of Flagging_Summaries directory
Flagged_fig, Flagged_plot = FlagStat_plot(FLAG_STAT, n)
thumb_img = 'thumb_'+ str(sizeX) + '_'+ Flagged_plot
make_Thumbnail(Flagged_fig, thumb_img, sizeX, sizeY, fig_dir)

# NOTE: careful - the location of Flagging_Summaries directory
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
			<th style="width:220"><img style="border:none;" src="" width="220" align="center"></th>
		    </tr></table>""")

html.write("""<h2 align="middle">Observation</h2>
                    <table width="100%" border="1">
                    <tr>
                        <th>SBID</th>
                        <th>No. of Antennas</th>
                        <th>Obs Start Date/Time<br>(UTC)</th> 
                        <th>Obs End Date/Time<br>(UTC)</th>   
                        <th>Duration<br>(hr)</th>
                        <th>Field</th>
                        <th>R.A.</th>
                        <th>Decl.</th>
                        <th>Total Bandwidth <br>(MHz)</th>
                    </tr>
		    <tr></tr>
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
                        <th>ASKAPsoft<br>version</th>
                        <th>Cal SBID</th>
                        <th>Frequency Range<br>(MHz)</th>
                        <th>Central Frequency<br>(MHz)</th>
                        <th>No. of Channels</th>
                        <th>Channel Width<br>(kHz)</th>
                        <th>Synthesised Beam<br>(arcsec &times; arcsec)</th>
                        <th>Beam Logs</th>
                        <th>Flagged Fractions</th>
                        <th>Flagged Antennas</th>
                        <th>Expected RMS</th>
                    </tr>
                    <tr align="middle">
                        <td>{0}</td>
                        <td>{1}</td>
                        <td>{2}</td>
                        <td>{3}</td>
                        <td>{4}</td>
                        <td>{5}</td>
                        <td>{6} &times; {7}</td>
                        <td>
                        <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                          <form action="{20}" method="get" target="_blank">
                          <button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
                         </form>
                        </td>
                        <td>
                        <a href="{21}" target="_blank"><img src="{22}" width="{23}" height="{24}" alt="thumbnail"></a>
                        """.format(askapsoft,
                                   cal_sbid,
                                   freq_range,                                   
                                   cfreq,
 				   nchan,
                                   chan_width_kHz,
                                   header_bmaj,
                                   header_bmin,
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

html.write("""</td>
                    </tr>
                    </table>
                    <br><table width="100%" border="1">
                    <tr>
                    <th colspan="4">Continuum Subtracted Beam Cube</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>MAD of Max Flux Density</p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>Mean RMS</p>
                    </td>
                    <td>
                    <a href="{8}" target="_blank"><img src="{9}" width="{10}" height="{11}" alt="thumbnail"></a>
                    <br><p>Median RMS</p>
                    </td>
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
                    <br><p>1-percentile noise rank</p>
                    """.format(fig_dir+'/'+'beamStat_MADMFD.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_MADMFD.png',
                               sizeX,
                               sizeY,
                               fig_dir+'/'+'beamStat_MeanRMS.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_MeanRMS.png',
                               sizeX,
                               sizeY,
                               fig_dir+'/'+'beamStat_MedianRMS.png',
                               fig_dir+'/'+ 'thumb_'+ str(sizeX) + '_beamStat_MedianRMS.png',
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
                    </td>
                    </tr>
                    <tr>
                    <th>Number of Bad Quality Channel</th>
                    <th>Number of Missing Data Channel</th>
                    <th>Bandpass Smoothing Issue*</th>
                    </tr>
                    <tr align="middle">
                    <td><br>{12}                  
                    <br><br><form action="{13}" method="get" target="_blank">
			<button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
                    </form>
		    </td>
                    <td><br>{14}
                    <br><br><form action="{15}" method="get" target="_blank">
			<button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
                    </form>
                    </td>
                    <td>{16}
                    </td>
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
                               QC_badchan_keyword,
                               fig_dir+'/' + str(mosaic_bad_chan),
                               QC_mdata_chan_keyword,
                               fig_dir+'/' + str(mosaic_mdata_chan),
                               bs_issue_id))


html.write("""</td> 
              </tr>
              </table>
              <p align=left>
              * SBIDs processed prior to March 2021 were affected by issues with bandpass smoothing parameterisation, resulting in a higher incidence of glitches on 1 MHz intervals.</p>
              <h2 align="middle">NVSS sources within 6 &times; 6 sq degree (Condon+ 1998)</h2>
              <table>
              <tr align="middle">
              <td id="cell">
              <form action="{0}" method="get" target="_blank">
                 <button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
              </form>
              """.format(fig_dir+'/' + nvss_cat))


##### Finish HTML report with generated time stamp
html.write("""
              </td>
              </tr>
              </table>                                             
              <p align=left>               
              Generated at {0} <br>
              <i> Report bugs to 
              <a href="mailto:hyein.yoon@sydney.edu.au">Hyein Yoon</a></i>
              </p>
              </body>
              </html>""".format(str(datetime.now())))


html.close()
print ("Spectral line validation report written to '{0}'.".format(html_name))
