################################################################################
# ASKAP-NAASHI data validation script
################################################################################

################################################################################
#
# This script generates an ASKAP-NAASHI spectral line cube HTML report.
# Compatibility: Python version 3.x
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- metadata/mslist-cal*txt
# -- metadata/mslist-*101.txt
# -- slurmOutput/*sh
# -- diagnostics/cubestats-*/*txt (any cubestats-* directory, e.g. cubestats-NAASHI-143)
# -- diagnostics/*png
# -- diagnostics/Flagging_Summaries/*averaged.ms.flagSummary
# -- SpectralCube_BeamLogs/*.txt
# -- image.restored.i.SB<SBID>.cube.contsub.fits (to read bmaj, bmin from header)
#
# Output files: all saved in Figures directory.
#
# To run:
#
# python <script name> -s <SBID> -c <cal_SBID>
# e.g. python askapnaashi_validation.py -s <SBID> -c <cal_SBID>
#      python askapnaashi_validation.py -s <SBID> -c <cal_SBID> -d /path/to/<SBID>
#
# ADMIN: all server paths / file names are in the ">>> ADMIN: EDIT PATHS" block below.
#
# Reference:
#
# 1) ASKAP-WALLABY data validation by Bi-Qing For
# Version: 24 March 2020
#
# 2) ASKAP-FLASH data validation by Hyein Yoon
# Original Version: 17 Nov 2021
# Adopted Version: 28 Feb 2024 (Github: askap-qc/validation/FLASH)
#
# 3) ASKAP-FLASH spectral line plot by Elizabeth Mahony
# Last update: 2 November 2026

# Adapted for ASKAP-NAASHI data validation (based on the FLASH version above)
# Last update: 24 November 2026
# Contact: Hyein Yoon
#
# NOTE: survey-specific values are collected in the "Survey configuration"
#       block below (confirmed for NAASHI; same band and setup as FLASH).
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

from scipy.constants import k as k_B
from scipy.stats import iqr
from scipy.optimize import curve_fit

from argparse import ArgumentParser

import matplotlib as mpl
mpl.use('Agg')  # to avoid matplotlib using the Xwindows backend

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.patches as mpatches


################################################################################
# >>> ADMIN: EDIT PATHS <<<
# Every server-specific path / file name is in this block (search for "ADMIN").
# Relative paths are relative to the SBID directory (DATA_DIR below).
################################################################################

# >>> ADMIN: EDIT PATH - SBID data directory. {sbid} is replaced by the SBID.
#     '.' = run the script from inside the SBID directory (old behaviour).
#     e.g. '/askapbuffer/naashi/{sbid}'. Can also be set with -d/--data_dir.
DATA_DIR = '.'

# >>> ADMIN: EDIT PATH - sub-directories inside the SBID directory
METADATA_DIR      = 'metadata'                 # mslist-*.txt, footprintOutput-*.txt, fieldlist*.txt
SLURM_DIR         = 'slurmOutput'              # *.sh (ASKAPsoft version)
DIAGNOSTICS_DIR   = 'diagnostics'              # cubestats-*/, Flagging_Summaries/
CUBESTATS_DIR_GLOB = 'cubestats-*'             # inside DIAGNOSTICS_DIR (e.g. cubestats-NAASHI_143)
FLAG_SUMMARY_GLOB = 'Flagging_Summaries/*averaged.ms.flagSummary'   # inside DIAGNOSTICS_DIR
BEAMLOG_DIR       = 'SpectralCube_BeamLogs'    # beamlog.image.restored.*.txt
SOURCE_SPECTRA_DIR = 'SourceSpectra'           # spec_SB{sbid}_component_{cid}.fits
NOISE_SPECTRA_DIR  = 'NoiseSpectra'            # noise_SB{sbid}_component_{cid}.fits

# >>> ADMIN: EDIT PATH - file name patterns ({sbid} -> SBID, {cid} -> component id)
MSLIST_GLOB         = 'mslist-*txt'                 # inside METADATA_DIR
MSLIST_SCIENCE_GLOB = 'mslist-scienceData*txt'      # inside METADATA_DIR
SLURM_GLOB          = '*.sh'                        # inside SLURM_DIR
CONTSUB_CUBE_FILE   = 'image.restored.{imagebase}.contsub.fits'   # optional (beam size from header)
SOURCE_SPECTRA_FILE = 'spec_SB{sbid}_component_{cid}.fits'
NOISE_SPECTRA_FILE  = 'noise_SB{sbid}_component_{cid}.fits'
SELAVY_CATALOGUE_GLOBS = ['casda/catalogues/selavy-image.i*.SB{sbid}.cont.*taylor.0.restored*.components.xml',
                          'catalogues/selavy-image.i*.SB{sbid}.cont.*taylor.0.restored*.components.xml',
                          'selavy*/selavy-image.i*.SB{sbid}.cont.*taylor.0.restored*.components.xml',
                          'selavy-image.i*.SB{sbid}.cont.*taylor.0.restored*.components.xml']   # first match is used
CONT_IMAGE_GLOBS = ['image.i.*SB{sbid}.cont.taylor.0.restored.conv.fits',
                    'image.i.*SB{sbid}.cont.taylor.0.restored*.fits']
FOOTPRINT_GLOBS = [METADATA_DIR + '/footprintOutput-sb{sbid}-src1-*.txt',
                   METADATA_DIR + '/footprintOutput-sb{sbid}-*.txt']
FIELDLIST_GLOB = METADATA_DIR + '/fieldlist*.txt'   # field centre; falls back to the mslist position

# >>> ADMIN: EDIT PATH - outputs (relative to DATA_DIR, or absolute)
FIG_DIR   = 'Figures'       # all plots / txt files
HTML_NAME = 'index.html'    # the report


################################################################################
# Survey configuration (NAASHI)
################################################################################

SURVEY_NAME = 'NAASHI'
REPORT_TITLE = 'ASKAP NAASHI Spectral Line Data Validation Report'
CONTACT_NAME = 'Hyein Yoon'
CONTACT_EMAIL = 'hiyoon.astro@gmail.com'
N_BLOCKS = 6          # correlator cards used to compute total bandwidth (FLASH pilot value)

# Theoretical RMS parameters
TSYS = 60             # K  Band 1 (NAASHI uses the same 711.5-999.5 MHz band as FLASH)
ANT_DIAM = 12         # metre
APER_EFF = 0.7        # aperture efficiency
COR_EFF = 0.8         # correlator efficiency
N_POL = 2.0           # number of polarisations

# QC thresholds
BEAMLOG_TOLERANCE = 0.3     # fractional deviation from median bmaj/bmin
BEAMLOG_MAX_OUTLIERS = 50   # channels allowed outside tolerance before a beam fails
BADCHAN_NSIGMA = 3.0        # bad channel threshold (x std of channel rms)
NOISERANK_SIGMA = 0.5       # Gaussian sigma for 1-percentile noise rank

# SBIDs affected by the bandpass smoothing issue (FLASH-specific list, prior to Mar 2021).
# Empty for NAASHI.
BS_ISSUE_SBIDS = []

# Catalogue Statistics cell (Mosaic Statistics table)
CAT_BRIGHT_MJY = 30.         # count Selavy components with flux_int above this (mJy); the pipeline extracts spectra by integrated flux

# Reference catalogue for the source list in the report (queried from VizieR)
#   'RACS-low' : RACS-low DR1, 887.5 MHz (Hale+ 2021), Dec < ~+30 deg
#   'RACS-mid' : RACS-mid DR1, 1367.5 MHz (Duchesne+ 2024)
#   'NVSS'     : NVSS, 1.4 GHz (Condon+ 1998), Dec > -40 deg
REF_CATALOGUE = 'RACS-low'
REF_CATALOGUES = {
    'RACS-low': {'vizier': 'J/other/PASA/38.58/gausscut', 'label': 'RACS-low',
                 'ref': 'Hale+ 2021', 'outfile': 'racs_low.txt'},
    'RACS-mid': {'vizier': 'J/other/PASA/41.3', 'label': 'RACS-mid',
                 'ref': 'Duchesne+ 2024', 'outfile': 'racs_mid.txt'},
    'NVSS':     {'vizier': 'VIII/65/nvss', 'label': 'NVSS',
                 'ref': 'Condon+ 1998', 'outfile': 'nvss.txt'},
}
# Source spectra (shown below Mosaic Statistics)
# Plotting follows plot_spec_local.py (E. Mahony): flux spectrum in 8 stacked panels,
# +/-5 sigma noise band from the noise spectrum, beam-forming intervals in blue,
# title 'Component <id> - <name> (<peak flux> Jy)' from the Selavy component catalogue.
# Components 1a-15a are plotted (Selavy numbers islands in order of decreasing flux).
# Input paths: see the ADMIN: EDIT PATHS block above
SOURCE_COMPONENTS = ['%da' % i for i in range(1, 16)]   # 1a, 2a, ..., 15a (3 rows of 5)
N_SOURCE_SPECTRA = len(SOURCE_COMPONENTS)
# Beam-forming intervals (blue lines): FLASH values, 711.5 MHz edge + every 5 MHz from 713.5 MHz
BF_INTERVAL_MHZ = 5.0       # NAASHI band = FLASH band (711.5-999.5 MHz)
BF_FIRST_MHZ = 713.5
SPECTRA_NCOLS = 5           # 5 columns per row
SPECTRA_THUMB_SIZE = 250    # thumbnail size in pixels
N_SPEC_PANELS = 8           # number of stacked frequency panels in each spectrum figure
SPEC_NOISE_NSIGMA = 5.0     # half-width of the grey noise band (x noise), as in plot_spec_local.py
SPEC_YLIM_NSIGMA = 20.0     # y-range = +/- this x median noise, as in plot_spec_local.py
SPEC_DPI = 150              # plot_spec_local.py uses 300

# Continuum image + Selavy bright components + beam footprint overlay (shown below Mosaic Statistics)
# Adapted from 06_askap_first_overlays_footprint.py (H. Yoon)
# Input paths: see the ADMIN: EDIT PATHS block above
CONT_VMIN, CONT_VMAX = 0.0, 0.003     # Jy/beam colour scale
CONT_IMAGE_STEP = 4                   # use every Nth pixel for plotting (speed/memory). 1 = full resolution
OVERLAY_FLUX_CUTS = (100., 200., 300.)  # mJy/beam: black 100-200, blue 200-300, red > 300

# Noise spectra statistics (adapted from scripts 01-05, H. Yoon)
INNER_RADIUS_DEG = 3.2         # inside / outside split (distance from the field centre)
DISTANCE_FLAT_SKY = False      # True: distance/offsets as in the original scripts 01-05,
                               #       sqrt(dRA^2 + dDec^2) without cos(Dec) (to compare with old results)
NOISE_YLIM = (0.001, 0.012)    # Jy/beam, rms vs channel plots
NOISE_DIST_YLIM = (0.001, 0.015)
NOISE_SKY_VLIM = (0.003, 0.01) # Jy/beam, colour scale of the spectral noise sky plot
CONT_NOISE_SKY_VLIM = (0.0, 1.0)  # mJy/beam, colour scale of the continuum noise sky plot
SKY_PLOT_HALFWIDTH = 4.5       # deg
NOISE_PLOT_ALL = False         # True: save a noise-vs-channel plot for every component in Figures/noise (slow, not shown in report)
NOISE_PLOT_RC = {'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 12,
                 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11}
NOISE_DIR = 'noise'            # sub-directory of Figures for the noise outputs
CONT_NOISE_COLUMNS = ['rms_image']   # Selavy continuum noise column (mJy/beam)

QUERY_CATALOGUE = False     # RACS query is slow - off by default. Set True (or run with --catalogue) to query VizieR
# NOTE: the RACS-low section is currently commented out in the main program (lines marked '# [RACS] ').

# Reference catalogue field width in degrees (box of width x width)
CAT_FIELD_WIDTH_DEG = 6


################################################################################
# Functions for the main program
################################################################################

#-------------------------------------------------------------------------------
# Find input files (stop with a clear message if missing)
#-------------------------------------------------------------------------------
def find_files(pattern, min_n=1):
    """
    Sorted glob that stops with a clear message if too few files are found.
    """

    files = sorted(glob.glob(pattern))
    if len(files) < min_n:
        sys.exit("ERROR: expected at least %d file(s) matching '%s' in %s, found %d.\n"
                 "       Run this script from inside the SBID directory."
                 % (min_n, pattern, os.getcwd(), len(files)))
    return files


#-------------------------------------------------------------------------------
# ASKAP 36-beam positions for the beam-layout plots
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# Gaussian model (for the noise-rank fit)
#-------------------------------------------------------------------------------
def gauss(x, *p):
    """
    Fitting a Gaussian function.
    """

    A, mu, sigma = p

    return A*np.exp(-(x-mu)**2.0/(2.0*sigma**2.0))


#-------------------------------------------------------------------------------
# Read beam size (bmaj, bmin) from the cube FITS header
#-------------------------------------------------------------------------------
def get_FitsHeader(fitsimage):
    """
    Getting basic information of processed mosaic.
    """

    hdu = getheader(fitsimage)
    bmaj = round(float(hdu['BMAJ'])*3600., 2)  # arcsec
    bmin = round(float(hdu['BMIN'])*3600., 2)  # arcsec

    return bmaj, bmin


#-------------------------------------------------------------------------------
# Convert sexagesimal RA/Dec strings to degrees
#-------------------------------------------------------------------------------
def radec_to_deg(ra, dec):
    """
    Converting mslist RA (hh:mm:ss.s) and Dec (dd.mm.ss.s) strings to decimal degrees.
    """

    TOKS_RA = ra.split(":")
    ra_deg = 15.0*(float(TOKS_RA[0]) + float(TOKS_RA[1])/60. + float(TOKS_RA[2])/3600.)

    sign = -1.0 if dec.strip().startswith('-') else 1.0
    TOKS_DEC = dec.strip().lstrip('+-').split(".", 2)
    dec_deg = sign*(float(TOKS_DEC[0]) + float(TOKS_DEC[1])/60. + float(TOKS_DEC[2])/3600.)

    return round(ra_deg, 5), round(dec_deg, 5)


#-------------------------------------------------------------------------------
# Reference catalogue query (RACS-low, VizieR) - currently disabled
#-------------------------------------------------------------------------------
def get_RefCatalogue(ra, dec):
    """
    Getting reference catalogue sources (default: RACS-low) within the
    CAT_FIELD_WIDTH_DEG x CAT_FIELD_WIDTH_DEG field through VizieR.
    The catalogue is chosen with REF_CATALOGUE in the Survey configuration.
    """

    from astroquery.vizier import Vizier # for installation: conda install -c astropy astroquery

    cat = REF_CATALOGUES[REF_CATALOGUE]
    catalogue = cat['vizier']
    out_file = cat['outfile']

    print ("Retrieving " + cat['label'] + " sources from Vizier. Depending on server connection, this might take a while....")

    ra_deg, dec_deg = radec_to_deg(ra, dec)
    centre = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='icrs')

    # Try the table ID first, then the parent catalogue, then all columns.
    # VizieR can return None (astroquery) instead of an empty list when a query is not understood.
    parent = '/'.join(catalogue.split('/')[:4]) if catalogue.startswith('J/') else catalogue
    attempts = [(catalogue, ['*', '_RAJ2000', '_DEJ2000']), (catalogue, ['**']),
                (parent, ['*', '_RAJ2000', '_DEJ2000'])]
    result, errors = None, []
    for cat_id, cols in attempts:
        try:
            v = Vizier(columns=cols, catalog=cat_id, timeout=600)
            v.ROW_LIMIT = -1
            result = v.query_region(centre, width=CAT_FIELD_WIDTH_DEG*u.deg)
        except Exception as e:
            errors.append('%s %s: %s' % (cat_id, cols, e))
            result = None
            continue
        if result is not None and len(result) > 0:
            print ("VizieR query OK: %s, columns %s, tables %s" % (cat_id, cols, list(result.keys())))
            break
        errors.append('%s %s: %s' % (cat_id, cols, 'no result' if result is None else 'no rows'))
        result = None

    if result is None:
        msg = "No %s sources found / query failed:\n  " % cat['label'] + "\n  ".join(errors)
        print (msg)
        print (msg, file=open(fig_dir + '/' + out_file, 'w'))
        return out_file

    # largest table = component catalogue
    table = max(result, key=len)
    print ("Found %d %s sources" % (len(table), cat['label']))
    table.write(fig_dir + '/' + out_file, format='ascii.fixed_width', delimiter=' ', overwrite=True)

    return out_file

    if len(result) > 0:
        print ("Found %d %s sources" % (len(result[0]), cat['label']))
        result[0].write(fig_dir + '/' + out_file, format='ascii.fixed_width', delimiter=' ', overwrite=True)
    else:
        print ("No sources found", file=open(fig_dir + '/' + out_file, 'w'))

    return out_file


#-------------------------------------------------------------------------------
# ASKAPsoft version from the slurm .sh file
#-------------------------------------------------------------------------------
def get_Version(param):
    """
    Getting the latest ASKAPsoft version that is used for the data reduction.
    """

    line = subprocess.check_output(['grep', 'Processed with ASKAPsoft', param])
    str_line = line.decode('utf-8')

    newline = str_line.splitlines()[-1] #Get the most recent, in case there is more than one instance
    askapsoft = newline.split()[-1]

    return askapsoft


#-------------------------------------------------------------------------------
# Flagging summary: number of records, channels, expected count
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# Flagging summary: flagged fraction and fully flagged antennas per beam
#-------------------------------------------------------------------------------
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

    # ndmin=2 keeps this working when there are 0 or 1 fully flagged baselines
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)   # empty file = no fully flagged baselines
        data = np.loadtxt('temp.dat', ndmin=2)
    if data.size == 0:
        data = np.zeros((0, 3))
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

    if total_flagged_ant >= 1:
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


#-------------------------------------------------------------------------------
# Observation metadata from the mslist file
#-------------------------------------------------------------------------------
def get_Metadata(metafile):
    """
    Getting basic information on observed field (one field only).
    """

    mslist_file = open(metafile, 'r')
    LINES = mslist_file.readlines()
    mslist_file.close()

    nBlocks = N_BLOCKS  # number of correlator cards (see Survey configuration)

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

    # Known bandpass smoothing issue SBIDs (see Survey configuration)
    bs_issue_sbid = BS_ISSUE_SBIDS

    bs_issue_id="No"
    for j in range(len(bs_issue_sbid)):
        if bs_issue_sbid[j] == float(sbid):
            bs_issue_id="Yes"
            break

    return n_ant, start_obs_date, end_obs_date, tobs, field, ra, dec, total_obs_bw, bs_issue_id


#-------------------------------------------------------------------------------
# Channel width, centre frequency and number of channels
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# Frequency range of the mosaic cube (cubeStats)
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# Sort key for component ids (1a, 1b, 2a, ...)
#-------------------------------------------------------------------------------
def component_sort_key(fname):
    """
    Sort key from Selavy component id in the file name, e.g. component_12b -> (12, 'b').
    Files without a component id go to the end.
    """

    import re
    m = re.search(r'component_(\d+)([a-z]*)', os.path.basename(fname))
    if m:
        return (int(m.group(1)), m.group(2), fname)
    return (10**9, '', fname)


#-------------------------------------------------------------------------------
# Read a 1D spectrum from a FITS file
#-------------------------------------------------------------------------------
def read_Spectrum(fitsfile):
    """
    Reading a 1D spectrum FITS file. Returns frequency (MHz), flux and header.
    """

    hdu = fits.open(fitsfile)
    hdr = hdu[0].header
    flux = np.squeeze(hdu[0].data).astype(float)
    hdu.close()

    faxis = None
    for k in range(1, hdr['NAXIS']+1):
        if str(hdr.get('CTYPE%d' % k, '')).startswith('FREQ'):
            faxis = k
            break
    chan = np.arange(flux.size)
    if faxis is None:
        raise ValueError('no FREQ axis in ' + fitsfile)
    freq = (hdr['CRVAL%d' % faxis] + (chan + 1 - hdr['CRPIX%d' % faxis])*hdr['CDELT%d' % faxis])/1e6  # MHz

    return freq, flux, hdr


#-------------------------------------------------------------------------------
# Running rms of a spectrum
#-------------------------------------------------------------------------------
def running_rms(flux, width):
    """
    Running robust rms (1.4826 x MAD). Only used when no noise spectrum exists.
    """

    rms = np.full(flux.size, np.nan)
    half = width//2
    for i in range(0, flux.size, half):
        seg = flux[max(0, i-half):i+half]
        seg = seg[np.isfinite(seg)]
        if seg.size > 10:
            rms[i:i+half] = 1.4826*np.median(np.abs(seg - np.median(seg)))
    return rms


#-------------------------------------------------------------------------------
# Selavy catalogue: component name and peak flux (for spectra titles)
#-------------------------------------------------------------------------------
def read_Selavy(sbid):
    """
    Reading component name and peak flux (Jy) from the Selavy component catalogue.
    Returns {component number (e.g. '1a'): (name, peak flux in Jy)}.
    """

    from astropy.io import votable

    cat_dict = {}
    for pattern in SELAVY_CATALOGUE_GLOBS:
        files = sorted(glob.glob(pattern.format(sbid=sbid)))
        if files:
            table = votable.parse_single_table(files[0]).array
            for line in table:
                compno = str(line['component_id']).split('_')[-1]
                compname = line['component_name']
                if isinstance(compname, bytes):
                    compname = compname.decode()
                cat_dict[compno] = (str(compname).strip(), float(line['flux_peak'])/1000.)  # mJy -> Jy
            print('Read %d components from %s' % (len(cat_dict), files[0]))
            break
    else:
        print('No Selavy component catalogue found - spectrum titles will have no peak flux')

    return cat_dict


#-------------------------------------------------------------------------------
# Beam footprint (beam centres) from footprintOutput
#-------------------------------------------------------------------------------
def read_Footprint(footprint_file):
    """
    Reading beam centres from the footprintOutput file
    (lines like '0 (-0.000 0.000) 12:30:00.000,-45:00:00.00').
    """

    beam_ra, beam_dec, beam_num = [], [], []
    with open(footprint_file) as f:
        lines = [l for l in f if ')' in l]
    for j, line in enumerate(lines):
        radec = line.split(')')[1].split(',')
        c = coord.SkyCoord(radec[0].strip(), radec[1].strip(), unit=(u.hour, u.deg), frame='icrs')
        beam_ra.append(c.ra.deg)
        beam_dec.append(c.dec.deg)
        beam_num.append(str(j))

    return np.array(beam_ra), np.array(beam_dec), beam_num


#-------------------------------------------------------------------------------
# Plot: continuum image + bright Selavy components + beam footprint
#-------------------------------------------------------------------------------
def ContinuumOverlay_plot(sbid, field, tobs_hr):
    """
    Continuum image with Selavy bright components and the beam footprint,
    adapted from 06_askap_first_overlays_footprint.py (H. Yoon).
    """

    from astropy.wcs import WCS
    from astropy.io import votable

    def first(patterns):
        for pattern in patterns:
            files = sorted(glob.glob(pattern.format(sbid=sbid)))
            if files:
                return files[0]
        return None

    image_file = first(CONT_IMAGE_GLOBS)
    if image_file is None:
        print('No continuum image found - skipping continuum overlay')
        return None, None
    component_file = first(SELAVY_CATALOGUE_GLOBS)
    footprint_file = first(FOOTPRINT_GLOBS)

    # continuum image (memory-mapped, every CONT_IMAGE_STEP-th pixel)
    hdul = fits.open(image_file, memmap=True)
    hdr = hdul[0].header
    data = hdul[0].data
    while data.ndim > 2:
        data = data[0]
    step = max(1, int(CONT_IMAGE_STEP))
    img = np.array(data[::step, ::step], dtype=np.float32)
    hdul.close()
    wcs = WCS(hdr, naxis=2).slice((slice(None, None, step), slice(None, None, step)))

    fig = plt.figure(figsize=(9, 8))
    fig.subplots_adjust(top=0.92, bottom=0.08, left=0.03, right=0.92)
    ax = fig.add_subplot(111, projection=wcs)
    im = ax.imshow(img, origin='lower', cmap='viridis_r', vmin=CONT_VMIN, vmax=CONT_VMAX)

    # beam footprint
    if footprint_file:
        bra, bdec, bnum = read_Footprint(footprint_file)
        ax.scatter(bra, bdec, transform=ax.get_transform('fk5'), s=20000, facecolors='none',
                   edgecolors='gray', lw=0.5, zorder=1)
        for i in range(len(bnum)):
            ax.text(bra[i]+0.08, bdec[i], bnum[i], transform=ax.get_transform('fk5'), fontsize=12, color='gray')
    else:
        print('No footprint file found - beams not drawn')

    # Selavy bright components
    if component_file:
        cat = votable.parse_single_table(component_file).array
        ra = np.array(cat['ra_deg_cont'], dtype=float)
        dec = np.array(cat['dec_deg_cont'], dtype=float)
        fpeak = np.array(cat['flux_peak'], dtype=float)   # mJy/beam
        compno = np.array([str(c.decode() if isinstance(c, bytes) else c).split('_')[-1] for c in cat['component_id']])

        c1, c2, c3 = OVERLAY_FLUX_CUTS
        for mask, size, color, lw in [((fpeak > c3), 50, 'red', 2),
                                      ((fpeak > c2) & (fpeak < c3), 35, 'blue', 1),
                                      ((fpeak > c1) & (fpeak < c2), 20, 'black', 1)]:
            ax.scatter(ra[mask], dec[mask], transform=ax.get_transform('fk5'), s=size,
                       facecolors='none', edgecolors=color, lw=lw, zorder=10)
            for o in np.where(mask)[0]:
                ax.text(ra[o]-0.05, dec[o]+0.05, compno[o], transform=ax.get_transform('fk5'), fontsize=7, color=color)

        ax.text(0.04, 0.02, 'f$_{peak}$ > %d mJy/beam' % c3, color='red', fontsize=10, transform=ax.transAxes)
        ax.text(0.34, 0.02, '%d < f$_{peak}$ < %d mJy/beam' % (c2, c3), color='blue', fontsize=10, transform=ax.transAxes)
        ax.text(0.68, 0.02, '%d < f$_{peak}$ < %d mJy/beam' % (c1, c2), color='black', fontsize=10, transform=ax.transAxes)
    else:
        print('No Selavy component catalogue found - components not drawn')

    ax.set_title('SB%s (%s, %.0fhr): Distribution of Bright Components (Continuum Image)' % (sbid, field, tobs_hr), fontsize=12)
    ax.grid(True, linestyle='--')
    ax.set_xlim(0, img.shape[1])
    ax.set_ylim(0, img.shape[0])
    ax.set_xlabel('R.A. (J2000)')
    ax.set_ylabel('Dec. (J2000)')

    cbar_ax = fig.add_axes([0.86, 0.08, 0.02, 0.84])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label('Jy/beam')

    plot_name = 'footprint_continuum_%s_SB%s.png' % (field, sbid)
    saved_fig = fig_dir + '/' + plot_name
    fig.savefig(saved_fig, dpi=150)
    plt.close(fig)

    return saved_fig, plot_name


#-------------------------------------------------------------------------------
# Field centre from fieldlist (fallback: mslist position)
#-------------------------------------------------------------------------------
def get_FieldCentre(ra, dec):
    """
    Field centre (deg) from metadata/fieldlist*.txt (last line: <id> <name> hh:mm:ss dd.mm.ss.s),
    falling back to the mslist position.
    """

    files = sorted(glob.glob(FIELDLIST_GLOB))
    if files:
        lines = [l.split() for l in open(files[0]) if l.strip()]
        try:
            return radec_to_deg(lines[-1][2], lines[-1][3])
        except Exception as e:
            print('Could not read field centre from ' + files[0] + ': ' + str(e))
    return radec_to_deg(ra, dec)


#-------------------------------------------------------------------------------
# Angular distance and RA/Dec offsets from the field centre
#-------------------------------------------------------------------------------
def sky_offset(ra, dec, ra_c, dec_c):
    """
    RA/Dec offsets (deg) from the field centre, with RA wrap-around at 0/360 and cos(Dec),
    and the angular distance (deg).
    """

    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    ddec = dec - dec_c
    if DISTANCE_FLAT_SKY:
        dra = (ra - ra_c + 180.) % 360. - 180.
        return dra, ddec, np.sqrt(dra**2 + ddec**2)
    dra = ((ra - ra_c + 180.) % 360. - 180.)*np.cos(np.radians(dec_c))
    dist = coord.SkyCoord(ra*u.deg, dec*u.deg).separation(coord.SkyCoord(ra_c*u.deg, dec_c*u.deg)).deg
    return dra, ddec, dist


#-------------------------------------------------------------------------------
# Selavy catalogue: position, flux_peak, rms_image, flux_int per component
#-------------------------------------------------------------------------------
def read_Selavy_table(sbid):
    """
    Selavy component table as {component_id: (ra, dec, flux_peak mJy, rms_image mJy, flux_int mJy)}.
    """

    from astropy.io import votable

    for pattern in SELAVY_CATALOGUE_GLOBS:
        files = sorted(glob.glob(pattern.format(sbid=sbid)))
        if files:
            cat = votable.parse_single_table(files[0]).to_table(use_names_over_ids=True)
            names = cat.colnames
            print('Selavy catalogue for noise statistics: %s (%d components)' % (files[0], len(cat)))
            noisecol = next((c for c in CONT_NOISE_COLUMNS if c in names), None)
            if noisecol is None:
                print('WARNING: none of %s in the Selavy catalogue - continuum noise cannot be plotted. '
                      'Columns: %s' % (CONT_NOISE_COLUMNS, ', '.join(names)))
            else:
                print('Continuum noise column: ' + noisecol)
            out = {}
            for line in cat:
                cid = line['component_id']
                cid = cid.decode() if isinstance(cid, bytes) else str(cid)
                rmsimg = float(line[noisecol]) if noisecol else np.nan
                out[cid.strip()] = (float(line['ra_deg_cont']), float(line['dec_deg_cont']),
                                    float(line['flux_peak']), rmsimg,
                                    float(line['flux_int']) if 'flux_int' in names else np.nan)
            return out
    return {}


#-------------------------------------------------------------------------------
# Noise spectra statistics (script 01): medians per third, distances
#-------------------------------------------------------------------------------
def NoiseSpectra_stats(sbid, ra_c, dec_c, out_dir):
    """
    Script 01: median noise of every noise spectrum (all channels and first/mid/last thirds),
    distance from the field centre, median text files, and noise-vs-channel plots.
    """

    files = sorted(glob.glob(os.path.join(NOISE_SPECTRA_DIR, NOISE_SPECTRA_FILE.format(sbid=sbid, cid='*'))),
                   key=component_sort_key)
    rows = []
    for name in files:
        hdu = fits.open(name)
        hdr = hdu[0].header
        noisedata = np.squeeze(hdu[0].data).astype(float)
        hdu.close()
        from astropy.wcs import WCS
        sky = WCS(hdr, naxis=2).pixel_to_world(0, 0)
        cid = component_sort_key(name)
        compno = '%d%s' % (cid[0], cid[1]) if cid[0] < 10**9 else os.path.basename(name)
        dist = sky_offset(sky.ra.deg, sky.dec.deg, ra_c, dec_c)[2]
        n3 = noisedata.size//3
        segs = [noisedata, noisedata[:n3], noisedata[n3:2*n3], noisedata[2*n3:3*n3]]
        meds = []
        for seg in segs:
            good = seg[np.isfinite(seg) & (seg < 1)]
            meds.append(np.median(good) if good.size else np.nan)
        rows.append(dict(compid='SB%s_component_%s' % (sbid, compno), compno=compno, file=name,
                         dist=float(dist), meds=meds, data=noisedata, n3=n3))

    if not rows:
        return None

    # median text files (component id, distance [deg], median noise [Jy/beam])
    txt_names = ['noise_median_all.txt', 'noise_median_first_third.txt',
                 'noise_median_mid_third.txt', 'noise_median_last_third.txt']
    for k, tname in enumerate(txt_names):
        with open(os.path.join(out_dir, tname), 'w') as f:
            f.write('# component_id  distance_from_field_centre(deg)  median_noise(Jy/beam)\n')
            for r in rows:
                f.write('%s %5.2f %5.9f\n' % (r['compid'], r['dist'], r['meds'][k]))

    # median over all components (per third), used as reference lines
    allmed = [np.nanmedian([r['meds'][k] for r in rows]) for k in range(4)]

    # noise-vs-channel plots
    plot_ids = None if NOISE_PLOT_ALL else set()   # per-component plots off unless NOISE_PLOT_ALL
    plots = {}
    for r in rows:
        if plot_ids is not None and r['compno'] not in plot_ids:
            continue
        n3 = r['n3']
        inside = r['dist'] < INNER_RADIUS_DEG
        color = 'blue' if inside else 'red'
        fig = plt.figure(figsize=(10, 3))
        fig.subplots_adjust(top=0.88, bottom=0.17, left=0.1, right=0.95)
        ax = fig.add_subplot(111)
        ax.plot(r['data'], '--', lw=0.7, color=color, alpha=0.5, zorder=-1)
        for k in range(3):
            ax.plot([k*n3, (k+1)*n3], [allmed[k+1]]*2, '--', color='orange', zorder=20)
            ax.text(k*n3 + 0.1*n3, 0.85, 'median= %.6f' % r['meds'][k+1], fontsize=9,
                    transform=ax.get_xaxis_transform())
        ax.set_title('%s (%s: d %s %.1f deg)' % (r['compid'], 'inside' if inside else 'outside',
                                                 '<' if inside else '>', INNER_RADIUS_DEG), color=color)
        ax.set_xlabel('channel')
        ax.set_ylabel('rms noise (Jy/beam)')
        ax.set_xlim(0, 3*n3)
        if allmed[1] > 0:
            ax.set_ylim(NOISE_YLIM[0], 2*allmed[1])
        pname = NOISE_DIR + '/noise_SB%s_component_%s.png' % (sbid, r['compno'])
        fig.savefig(fig_dir + '/' + pname, dpi=100)
        plt.close(fig)
        plots[r['compno']] = pname

    return dict(rows=rows, txt=txt_names, allmed=allmed, plots=plots, n3=rows[0]['n3'])


#-------------------------------------------------------------------------------
# Plot: spectral noise: inside vs outside (script 02)
#-------------------------------------------------------------------------------
def NoiseInOut_plot(sbid, stats):
    """
    Script 02: noise of all components per third, inside / outside / both.
    """

    pylab.rcParams.update(NOISE_PLOT_RC)

    rows, n3 = stats['rows'], stats['n3']
    dist = np.array([r['dist'] for r in rows])
    meds = np.array([r['meds'][1:] for r in rows])   # (ncomp, 3)
    inside = dist < INNER_RADIUS_DEG
    outside = ~inside

    def med(mask, k):
        v = meds[mask & (meds[:, k] > 0), k]
        return np.median(v) if v.size else np.nan

    fig = plt.figure(figsize=(9, 9))
    fig.subplots_adjust(top=0.95, bottom=0.07, left=0.13, right=0.95, hspace=0.65)
    panels = [(inside, 'blue', r'SB%s: all components - inside (d$_{\rm fieldcentre}$ < %.1f deg)' % (sbid, INNER_RADIUS_DEG)),
              (outside, 'red', r'all components - outside (d$_{\rm fieldcentre}$ > %.1f deg)' % INNER_RADIUS_DEG)]
    for p, (mask, color, title) in enumerate(panels):
        ax = fig.add_subplot(3, 1, p+1)
        for k in range(3):
            for v in meds[mask, k]:
                ax.plot([k*n3, (k+1)*n3], [v, v], '-', color='gold', lw=0.5, zorder=10)
            m = med(mask, k)
            ax.plot([k*n3, (k+1)*n3], [m, m], '--', color=color, zorder=20)
            ax.text(k*n3 + 0.5*n3, 0.88, 'median= %.8f' % m, fontsize=11, color=color, zorder=50,
                    ha='center', va='center', transform=ax.get_xaxis_transform())   # no background box; zorder keeps it above the lines
        ax.set_xlim(0, 3*n3)
        ax.set_ylim(*NOISE_YLIM)
        ax.set_title(title + '  (n=%d)' % mask.sum())
        ax.set_xlabel('channel')
        ax.set_ylabel('rms noise (Jy/beam)')
    ax = fig.add_subplot(3, 1, 3)
    for mask, color in [(inside, 'blue'), (outside, 'red')]:
        for k in range(3):
            m = med(mask, k)
            ax.plot([k*n3, (k+1)*n3], [m, m], '--', color=color, zorder=20)
    ax.set_xlim(0, 3*n3)
    ax.set_ylim(*NOISE_YLIM)
    ax.set_title('all components - both')
    ax.set_xlabel('channel')
    ax.set_ylabel('rms noise (Jy/beam)')

    pname = NOISE_DIR + '/SB%s_all_components.png' % sbid
    fig.savefig(fig_dir + '/' + pname, dpi=120)
    plt.close(fig)
    return pname


#-------------------------------------------------------------------------------
# Plot: rms noise vs distance from field centre (script 03)
#-------------------------------------------------------------------------------
def NoiseDist_plot(sbid, stats):
    """
    Script 03: median noise vs distance from the field centre.
    """

    pylab.rcParams.update(NOISE_PLOT_RC)

    rows, n3 = stats['rows'], stats['n3']
    dist = np.array([r['dist'] for r in rows])
    fig = plt.figure(figsize=(6, 4))
    fig.subplots_adjust(top=0.92, bottom=0.12, left=0.15, right=0.95)
    ax = fig.add_subplot(111)
    for k, (color, lab) in enumerate([('orange', 'first'), ('green', 'mid'), ('purple', 'last')]):
        ax.plot(dist, [r['meds'][k+1] for r in rows], 'o', ms=3, color=color, alpha=0.5, zorder=10,
                label='%s %d chan' % (lab, n3))
    ax.axvline(INNER_RADIUS_DEG, ls='--', lw=0.7, color='gray', alpha=0.5, zorder=-1)
    ax.set_xlim(0, max(5, np.nanmax(dist)*1.05))
    ax.set_ylim(*NOISE_DIST_YLIM)
    ax.set_title('SB%s' % sbid)
    ax.set_xlabel('distance from field centre (deg)')
    ax.set_ylabel('rms noise (Jy/beam)')
    ax.legend(numpoints=1)
    pname = NOISE_DIR + '/SB%s_rms_vs_dist_fieldcentre.png' % sbid
    fig.savefig(fig_dir + '/' + pname, dpi=120)
    plt.close(fig)
    return pname


#-------------------------------------------------------------------------------
# Plot: spectral / continuum noise sky distribution (scripts 04-05)
#-------------------------------------------------------------------------------
def NoiseSky_plot(sbid, stats, selavy_tab, ra_c, dec_c, footprint, continuum=False):
    """
    Scripts 04/05: sky distribution of the median noise (mid third of the spectral noise,
    or the continuum rms_image from Selavy), with the beam footprint and histograms.
    """

    pylab.rcParams.update(NOISE_PLOT_RC)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    rows = [r for r in stats['rows'] if r['compid'] in selavy_tab]
    if not rows:
        print('No noise spectra matched the Selavy catalogue - skipping noise sky plot')
        return None, ''
    ra = np.array([selavy_tab[r['compid']][0] for r in rows])
    dec = np.array([selavy_tab[r['compid']][1] for r in rows])
    noise = np.array([r['meds'][2] for r in rows])            # mid third
    cont = np.array([selavy_tab[r['compid']][3] for r in rows])
    dra, ddec, dist = sky_offset(ra, dec, ra_c, dec_c)
    good = noise > 0
    value = cont if continuum else noise
    if continuum:
        ok = np.isfinite(cont[good])
        print('Continuum noise: %d of %d matched components have continuum noise%s' % (
            ok.sum(), good.sum(), (', range %.4g-%.4g mJy/beam (colour scale %g-%g)' % (
                np.nanmin(cont[good]), np.nanmax(cont[good]), *CONT_NOISE_SKY_VLIM)) if ok.any() else ''))
    print('Noise sky plot: %d noise spectra matched to Selavy components (of %d)' % (len(rows), len(stats['rows'])))

    fig = plt.figure(figsize=(9, 8))
    fig.subplots_adjust(top=0.92, bottom=0.08, left=0.1, right=0.84)
    ax = fig.add_subplot(111)

    if footprint is not None:
        bdra, bddec, _ = sky_offset(footprint[0], footprint[1], ra_c, dec_c)
        ax.scatter(bdra, bddec, s=14000, facecolors='none', edgecolors='gray', lw=0.3, zorder=-10)
        for i, num in enumerate(footprint[2]):
            ax.text(bdra[i]+0.08, bddec[i], num, fontsize=12, color='gray')

    r0 = INNER_RADIUS_DEG
    ax.plot([r0, r0, -r0, -r0, r0], [-r0, r0, r0, -r0, -r0], ls='--', c='red', zorder=-10, lw=1.0, alpha=0.7)
    vmin, vmax = CONT_NOISE_SKY_VLIM if continuum else NOISE_SKY_VLIM
    im = ax.scatter(dra[good], ddec[good], s=40, c=value[good], edgecolor='gray', cmap='jet',
                    vmin=vmin, vmax=vmax, linewidth=1.0, zorder=10, alpha=0.7)

    hw = SKY_PLOT_HALFWIDTH
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes('top', 1.3, pad=0.04, sharex=ax)
    axHisty = divider.append_axes('right', 1.3, pad=0.04, sharey=ax)
    out = good & (dist > r0)
    for mask, fc, alpha in [(good, 'gray', 0.4), (out, 'red', 0.3)]:
        axHistx.hist(dra[mask], bins=18, facecolor=fc, edgecolor='gray', range=(-hw, hw), alpha=alpha)
        axHisty.hist(ddec[mask], bins=18, orientation='horizontal', facecolor=fc, edgecolor='gray',
                     range=(-hw, hw), alpha=alpha)
    plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)
    axHistx.grid(True, linestyle='--')
    axHisty.grid(True, linestyle='--')
    ax.grid(True, linestyle='--', lw=0.5)
    ax.set_xlim(hw, -hw)
    ax.set_ylim(-hw, hw)

    what = 'Continuum Image' if continuum else 'Mid %d Chan' % stats['n3']
    axHistx.set_title('SB%s: %s Noise Sky Distribution - %s' % (sbid, 'Continuum' if continuum else 'Spectral', what), fontsize=12)
    ax.set_xlabel('RA offset from the field centre (deg)', fontsize=12)
    ax.set_ylabel('Dec offset from the field centre (deg)', fontsize=12)
    cbar_ax = fig.add_axes([0.87, 0.08, 0.02, 0.84])
    fig.colorbar(im, cax=cbar_ax,
                 label='rms noise (mJy/beam)' if continuum else 'median noise flux density (Jy/beam)')

    pname = NOISE_DIR + ('/continuum_noise_sky_SB%s.png' if continuum else '/noise_sky_SB%s.png') % sbid
    fig.savefig(fig_dir + '/' + pname, dpi=120)
    plt.close(fig)

    ntot = int(good.sum())
    nout = int(out.sum())
    summary = '%d spectra, %d (%.1f%%) at d &gt; %.1f deg' % (ntot, nout, 100.*nout/max(ntot, 1), r0)

    # combined table (script 04/05 *_combined_all.txt)
    if not continuum:
        with open(fig_dir + '/' + NOISE_DIR + '/SB%s_combined_all.txt' % sbid, 'w') as f:
            f.write('# component_id ra_deg dec_deg flux_peak(mJy/beam) median_noise_mid(Jy/beam) '
                    'continuum_noise(mJy/beam) distance(deg)\n')
            for i, r in enumerate(rows):
                f.write('%s %5.6f %5.6f %5.3f %5.9f %5.9f %5.3f\n' % (r['compid'], ra[i], dec[i],
                        selavy_tab[r['compid']][2], noise[i], cont[i], dist[i]))

    return pname, summary


#-------------------------------------------------------------------------------
# Plot: source spectrum (plot_spec_local.py style, credit: Elizabeth Mahony)
#-------------------------------------------------------------------------------
def SourceSpectrum_plot(fitsfile, compno, cat_dict, noisefile):
    """
    Plotting a source flux spectrum, adapted from make_plot() in plot_spec_local.py
    (E. Mahony; plotting code from J. Allison's miriad_pipe.py).
    8 stacked panels, frequency decreasing to the right, HI redshift on the top axis,
    grey +/- 5 sigma noise band, blue lines at beam-forming intervals.
    """

    hi_rest = 1420.40575177   # MHz

    freq, flux, hdr = read_Spectrum(fitsfile)
    if noisefile and os.path.isfile(noisefile):
        noiseflux = read_Spectrum(noisefile)[1]
    else:
        print('No noise spectrum for component %s - using running rms of the spectrum' % compno)
        noiseflux = running_rms(flux, 200)

    if compno in cat_dict:
        compname, peak_flux = cat_dict[compno]
        mktitle = 'Component %s - %s (%.3f Jy)' % (compno, compname, peak_flux)
    else:
        compname = str(hdr.get('OBJECT', '')).strip()
        mktitle = 'Component %s - %s' % (compno, compname)

    # Reverse x-axis (high frequency first)
    freq = freq[::-1]
    flux = flux[::-1]
    noiseflux = noiseflux[::-1]
    z = hi_rest/freq - 1.0

    # Beam-forming intervals
    fmin, fmax = np.nanmin(freq), np.nanmax(freq)
    chunks = np.hstack([fmin, np.arange(BF_FIRST_MHZ, fmax, BF_INTERVAL_MHZ), fmax])

    ymax = np.nanmedian(noiseflux)*SPEC_YLIM_NSIGMA
    ymin = -ymax

    # A4 portrait plot
    with plt.rc_context({'font.family': 'serif', 'font.size': 14,
                         'xtick.labelsize': 10, 'ytick.labelsize': 10}):
        fig = plt.figure(figsize=(8.3, 11.7))
        fig_pad = 0.5
        nplots = N_SPEC_PANELS
        fig.subplots_adjust(wspace=0, hspace=fig_pad)

        for subplot in range(0, nplots):
            N = len(freq)/nplots
            ind_min = int(subplot*N)
            ind_max = int((subplot+1)*N-1)

            sub_x = np.nan_to_num(freq[ind_min:ind_max], nan=0.)
            sub_y = np.nan_to_num(flux[ind_min:ind_max], nan=0.)
            sub_n = np.nan_to_num(noiseflux[ind_min:ind_max], nan=1.)   # NaN (flagged) channels -> full grey stripe
            sub_z = np.nan_to_num(z[ind_min:ind_max], nan=0.)

            ax1 = fig.add_subplot(nplots, 1, subplot+1)
            if subplot == nplots//2:
                ylabh = ax1.set_ylabel('Flux density (Jy/beam)', fontsize=14)
                ylabh.set_position((ylabh.get_position()[1], 1.0 + 0.5*fig_pad))
            if subplot == nplots-1:
                ax1.set_xlabel('Frequency (MHz)')
            if subplot == 0:
                plt.suptitle(mktitle, x=0.5, y=0.95, horizontalalignment='center', fontsize=14)

            xmin, xmax = max(sub_x), min(sub_x)
            ax1.set_xlim(xmin, xmax)
            if np.isnan(ymax):
                continue
            ax1.set_ylim(ymin, ymax)

            ax1.plot(sub_x, sub_y, 'k-', linewidth=0.5, zorder=3)
            ax1.axhline(color='r', linestyle='-', linewidth=0.5, zorder=2)
            ax1.fill_between(sub_x, -SPEC_NOISE_NSIGMA*sub_n, SPEC_NOISE_NSIGMA*sub_n,
                             facecolor=[0.8, 0.8, 0.8], edgecolor='none', zorder=1)
            for i in chunks:   # beam-forming intervals
                ax1.plot([i, i], [ymin, ymax], color='b', linestyle='-', linewidth=0.5, zorder=2)

            # redshift labels on the top axis
            label_inc = 0.01
            label_min = label_inc*np.floor(min(sub_z)/label_inc)
            label_max = label_inc*np.ceil(max(sub_z)/label_inc)
            labels = np.arange(label_min, label_max, label_inc)
            positions = [hi_rest/(1+label) for label in labels]
            ax1.get_xaxis().tick_bottom()
            ax2 = fig.add_axes(ax1.get_position(), sharey=ax1, frameon=False)
            ax2.xaxis.tick_top()
            ax2.set_xticks(positions)
            ax2.set_xticklabels(['%.2f' % label for label in labels])
            ax2.set_xlim(xmin, xmax)
            ax2.yaxis.set_visible(False)
            if subplot == 0:
                ax1.text(0.5, 1.3, r'$z_\mathrm{HI}$', horizontalalignment='center', fontsize=12, transform=ax1.transAxes)

            ax1.minorticks_on()
            ax1.tick_params(bottom=True, left=True, top=False, right=False, length=6, width=1, which='major', direction='in')
            ax1.tick_params(bottom=True, left=True, top=False, right=False, length=3, width=1, which='minor', direction='in')
            ax2.tick_params(bottom=False, left=False, top=True, right=False, length=6, width=1, which='major', direction='in')
            ax2.tick_params(bottom=False, left=False, top=False, right=False, length=6, width=1, which='minor', direction='in')

        plot_name = 'SB%s_component_%s_flux.png' % (sbid, compno)
        saved_fig = fig_dir + '/' + plot_name
        fig.savefig(saved_fig, dpi=SPEC_DPI)
        plt.close(fig)

    label = 'Component ' + compno + ('<br>' + compname if compname else '')

    return saved_fig, plot_name, label


#-------------------------------------------------------------------------------
# Make a thumbnail image for the HTML report
#-------------------------------------------------------------------------------
def make_Thumbnail(image, thumb_img, sizeX, sizeY, dir):
    """
    Making thumbnail images.
    """

    size = sizeX, sizeY
    im = Image.open(image)
    im.thumbnail(size)
    im.save(dir+ '/' + thumb_img)

    return thumb_img


#-------------------------------------------------------------------------------
# QC: missing data channels (Max == 0) in the mosaic cube
#-------------------------------------------------------------------------------
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
        QC_mdata_chan_keyword = 'n= 0'   # 'none' placeholder is not a channel
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

    return  n_mdata_chan, mosaic_mdata_chan, QC_mdata_chan_id, QC_mdata_chan_keyword


#-------------------------------------------------------------------------------
# QC: bad quality channels in the mosaic cube
#-------------------------------------------------------------------------------
def qc_Bad_Chans(infile, mad_rms, med_rms):
    """
    Checking for bad data channels in the mosaic cube.
    Printing the list of bad data channels.
    """

    data = np.loadtxt(infile)
    data_rms = data[:,3]
    data_median_rms=np.median(data_rms)
    data_std_rms=np.std(data_rms)
    threshold = np.std(data_rms)*BADCHAN_NSIGMA  # threshold = N * standard deviation from med_rms

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
        QC_badchan_keyword = 'n= 0'   # 'none' placeholder is not a channel
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

    return n_bad_chan, mosaic_bad_chan, QC_badchan_id, QC_badchan_keyword


#-------------------------------------------------------------------------------
# QC: beam size (bmaj, bmin) variation from the beam logs
#-------------------------------------------------------------------------------
def qc_BeamLogs():
    """
     Threshold - BEAMLOG_TOLERANCE (default 30%) deviation from median bmaj, bmin of each beam.
     Beam sizes can change greatly across a wide band.
     Checking with some plots - e.g. bmaj or bmin vs chan is recommended.
     Note that bmaj and bmin for the first few channels are always zero due to barycentric correction.
    """

    file_dir = BEAMLOG_DIR
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
            tolerance_maj=[bmaj_median - bmaj_median*BEAMLOG_TOLERANCE, bmaj_median + bmaj_median*BEAMLOG_TOLERANCE]
            tolerance_min=[bmin_median - bmin_median*BEAMLOG_TOLERANCE, bmin_median + bmin_median*BEAMLOG_TOLERANCE]
            MEDIAN_BMAJ.append(bmaj_median)
            MEDIAN_BMIN.append(bmin_median)

        # check bmaj
            outliers_bmaj = (bmaj < tolerance_maj[0]) | (bmaj > tolerance_maj[1])

            if np.count_nonzero(outliers_bmaj) > BEAMLOG_MAX_OUTLIERS:
                qc_BMAJ_label = 'fail'
            else:
                qc_BMAJ_label = 'pass'

        # check bmin
            outliers_bmin = (bmin < tolerance_min[0]) | (bmin > tolerance_min[1])

            if np.count_nonzero(outliers_bmin) > BEAMLOG_MAX_OUTLIERS:
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


#-------------------------------------------------------------------------------
# Beam statistic: MAD of max flux density
#-------------------------------------------------------------------------------
def cal_beam_MADMFD(infile):
    """
    Calculating the MAD of max flux density of each beam.
    """

    data = np.loadtxt(infile)
    maxfdensity = data[:,8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 3)

    return mad_maxfdensity


#-------------------------------------------------------------------------------
# Beam statistic: mean rms
#-------------------------------------------------------------------------------
def cal_beam_MeanRMS(infile):
    """
    Calculating the average RMS of each beam.
    """

    data = np.loadtxt(infile)
    rms = data[:,3]
    mean_rms = round(np.mean(rms), 3)

    return mean_rms


#-------------------------------------------------------------------------------
# Beam statistic: median rms
#-------------------------------------------------------------------------------
def cal_beam_MedianRMS(infile):
    """
    Calculating the median RMS of each beam.
    """

    data = np.loadtxt(infile)
    rms = data[:,3]
    median_rms = round(np.median(rms), 3)

    return median_rms


#-------------------------------------------------------------------------------
# Mosaic statistic: median and MAD of channel rms
#-------------------------------------------------------------------------------
def cal_mosaic_Stats(infile):
    """
    Calculating MAD RMS and median RMS for the mosaic cube.
    """

    data = np.loadtxt(infile)
    rms = data[:,3]

    med_rms = np.median(rms)
    mad_rms = round(median_absolute_deviation(rms), 3)

    return mad_rms, med_rms


#-------------------------------------------------------------------------------
# Theoretical rms noise (radiometer equation)
#-------------------------------------------------------------------------------
def cal_Theoretical_RMS(n_ant, tobs, chan_width):
    """
    Calculating the theoretical rms noise for ASKAP.
    Assuming natural weighting and not taking into account fraction of flagged data.
    """

    tsys = TSYS          # K (see Survey configuration)
    antdiam = ANT_DIAM   # metre
    aper_eff = APER_EFF  # aperture efficiency
    coreff = COR_EFF     # correlator efficiency
    npol = N_POL         # Number of polarisation, npol = 2 for images in Stokes I, Q, U, or V

    anteff = np.pi*(antdiam/2.0)**2.0*aper_eff
    SEFD = 2.0*k_B*1e26*tsys/anteff
    rms_jy = SEFD/(coreff*np.sqrt(npol*n_ant*(n_ant-1.0)*chan_width*tobs))

    return rms_jy


#-------------------------------------------------------------------------------
# Expected rms per beam, corrected for flagging
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# All diagnostics/cubestats-* directories
#-------------------------------------------------------------------------------
def cubestats_dirs():
    """
    All diagnostics/cubestats-* directories, the one named after the field first
    (directory names can differ from the mslist field name, e.g. cubestats-NAASHI-143).
    """

    dirs = [d for d in sorted(glob.glob(os.path.join(diagnostics_dir, CUBESTATS_DIR_GLOB))) if os.path.isdir(d)]
    preferred = [diagnostics_dir + '/cubestats-' + field, diagnostics_dir + '/cubestats-' + field.replace('_', '-')]
    return [d for d in preferred if d in dirs] + [d for d in dirs if d not in preferred]


#-------------------------------------------------------------------------------
# Find files in any cubestats-* directory
#-------------------------------------------------------------------------------
def cubestats_glob(pattern):
    """
    Files matching pattern (e.g. '*.png') in any cubestats directory.
    """

    files = []
    for d in cubestats_dirs():
        files += sorted(glob.glob(d + '/' + pattern))
    return files


#-------------------------------------------------------------------------------
# Directory and file prefix of the per-beam cubeStats files
#-------------------------------------------------------------------------------
def cubestats_beam_paths():
    """
    Directory and file prefix of the per-beam cubeStats files.
    """

    for d in cubestats_dirs():
        for b in ['/cubeStats-image.restored.' + imagebase + '.' + field,
                  '/cubeStats-image.restored.' + 'i.' + field + '.SB' + sbid + '.cube.' + field]:
            if glob.glob(d + b + '*.txt'):
                return d, b
        hits = sorted(glob.glob(d + '/cubeStats-image.restored*.beam00.contsub.txt'))
        if hits:
            return d, hits[0][len(d):].split('.beam00')[0]
    print('WARNING: no per-beam cubeStats files found in ' + str(cubestats_dirs()))
    return diagnostics_dir + '/cubestats-' + field, '/cubeStats-image.restored.' + imagebase + '.' + field


#-------------------------------------------------------------------------------
# Plot: flagged visibility fraction per beam
#-------------------------------------------------------------------------------
def FlagStat_plot(FLAGSTAT, n):
    """
    Plotting and visualising flagging statistics of 36 beams.
    """

    file_dir, basename = cubestats_beam_paths()

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


#-------------------------------------------------------------------------------
# Plot: number of fully flagged antennas per beam
#-------------------------------------------------------------------------------
def FlagAnt_plot(N_FLAG_ANT, n):
    """
    Plotting and visualising number of flagged (completely) antennas beam-by-beam.
    """

    title = 'No. of 100% flagged antennas'
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


#-------------------------------------------------------------------------------
# Plot: expected rms per beam
#-------------------------------------------------------------------------------
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
    cb.set_label('mJy/beam')
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


#-------------------------------------------------------------------------------
# Plot: beam log QC per beam
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# QC: noise rank from the 1-percentile spread
#-------------------------------------------------------------------------------
def qc_NoiseRank(spread):
    """
    Evaluating the 1-percentile noise rank distribution.
    Sigma of the Gaussian function is used as a metric.
    """

    sigma = NOISERANK_SIGMA # Gaussian distribution (default value, sigma=0.5)

    if (spread <= 3*sigma):
        qc_label = 'good'
    elif (3*sigma < spread < 6*sigma):
        qc_label = 'ok'
    else:
        qc_label = 'uncertain'

    return qc_label


#-------------------------------------------------------------------------------
# Plot: 1-percentile noise rank histograms
#-------------------------------------------------------------------------------
def NoiseRank_histplot(nchan):
    """
    Plotting and visualising 1-percentile noise level histogram of 36 beams.
    """

    ID_LABEL = []
    plot_name = 'beam_1pctile_hist_SB'+ sbid + '.png'
    saved_fig = fig_dir + '/' + plot_name
    file_dir, basename = cubestats_beam_paths()

    params = {'axes.labelsize': 6,
              'axes.titlesize':6,
              'xtick.labelsize':5,
              'ytick.labelsize':5,
              'font.size':6}

    pylab.rcParams.update(params)

    fig, axs = plt.subplots(6,6)
    fig.subplots_adjust(hspace = .004, wspace=.004, bottom=0.05)
    fig.text(0.5, 0.007, '1-percentile noise level (mJy/beam)', ha='center')
    fig.text(0.01, 0.5, 'log N', ha='center', rotation='vertical')

    axs = axs.ravel()

    for i in range(0,36):
        infile = file_dir + basename +'.beam%02d.contsub.txt'%(i)
        if os.path.isfile(infile):
            data = np.loadtxt(infile)
            onepctile = data[:,6][data[:,3]>0]
            median_val = np.median(onepctile)
            # if statement is needed to rule out really bad data without having to do the Gaussian fitting
            #   add an IQR test to exclude the case of IQR=0 (ie all zeros)
            if (median_val > 1000.0) or (median_val < -1000.0) or not iqr(onepctile)>0.:
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
                print(f'binning with width={bin_width}, giving {n_bins} bins')

                hist, bins = np.histogram(onepctile, bins=n_bins, range=(xmin_val-3, xmax_val+3))
                with np.errstate(divide='ignore'):  # ignore division of zero
                    N = np.log10(hist)   # get log N for y-axis
                    N[N == -inf] = 0

                xcenter = (bins[:-1] + bins[1:]) / 2
                ymax_val = np.max(N)
                median_val_x = np.median(x)

                # Fitting a Gaussian and use spread (sigma) as a metric
                guess=[ymax_val, median_val_x, 5.0]
                print(f'Fitting Gaussian to {guess} with N={N} and xcenter={xcenter}')
                try:
                    coeff, var_matrix = curve_fit(gauss, xcenter, N, guess)
                    spread = round(np.abs(coeff[2]), 3)
                    ID_LABEL.append(qc_NoiseRank(spread))

                    axs[i].bar(xcenter, N)
                    axs[i].plot(xcenter,gauss(xcenter,*coeff),'r-',lw=1)
                    axs[i].set_xlim(xmin_val-3, xmax_val+3)
                    axs[i].set_ylim(0, ymax_val+3)
                    axs[i].title.set_text('Beam%02d' %(i))
                except:
                    print(f'Could not fit Gaussian for beam {i:02d}')
                    ID_LABEL.append('bad')
                    axs[i].set_xlim(-1, 1)
                    axs[i].set_ylim(0, 3)
                    axs[i].title.set_text('Beam%02d' %(i))

        else:
            ID_LABEL.append('bad')

    plt.tight_layout()
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name, ID_LABEL


#-------------------------------------------------------------------------------
# Plot: noise rank QC per beam
#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# Plot: beam statistic (MAD max flux, mean rms, median rms) per beam
#-------------------------------------------------------------------------------
def BeamStat_plot(item, n):
    """
    Plotting and visualising statistics of 36 beams.
    """

    file_dir, basename = cubestats_beam_paths()

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
    cb.set_label('mJy/beam')
    plt.savefig(saved_fig,bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


################################################################################
# Main program where it calls all the functions
################################################################################


# ignore astropy warnings
warnings.simplefilter('ignore', AstropyWarning)

parser = ArgumentParser(description='Run NAASHI validation and produce an HTML report')
parser.add_argument('-s','--sbid', dest='sbid',required='true',help='Science SBID',type=str)
parser.add_argument('-c','--cal_sbid', dest='cal_sbid',required='true',help='Calibrator SBID',type=str)
parser.add_argument('-i','--imagebase', dest='imagebase',default='i.SB%s.cube',help='Base string for images [default=%default]',type=str)
# [RACS] parser.add_argument('--catalogue', dest='catalogue', action='store_true',help='Query the reference catalogue (RACS) from VizieR (slow)')
parser.add_argument('-d','--data_dir', dest='data_dir', default=None, help='SBID data directory (default: DATA_DIR in the script)', type=str)
# [RACS] parser.add_argument('--no_catalogue', dest='no_catalogue', action='store_true',help='Skip the reference catalogue (RACS) query')
options = parser.parse_args()

sbid = options.sbid
# Work inside the SBID directory, so every relative path in the ADMIN block resolves there
data_dir = (options.data_dir or DATA_DIR).format(sbid=sbid)
if not os.path.isdir(data_dir):
    sys.exit("ERROR: data directory '%s' not found (set DATA_DIR or use -d)" % data_dir)
os.chdir(data_dir)
print('Data directory: ' + os.getcwd())

fig_dir = FIG_DIR
cal_sbid = options.cal_sbid

# beam number
n = [26,25,24,23,22,21,27,10,9,8,7,20,28,11,3,1,6,19,29,12,2,0,5,18,30,13,14,15,4,17,31,32,33,34,35,16]

html_name = HTML_NAME

diagnostics_dir = DIAGNOSTICS_DIR

imagebase = options.imagebase
imagebase = imagebase.replace('%s',sbid)

if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Required files

metafile = find_files(os.path.join(METADATA_DIR, MSLIST_GLOB))[0]
metafile_science = find_files(os.path.join(METADATA_DIR, MSLIST_SCIENCE_GLOB))[0]
param_file = find_files(os.path.join(SLURM_DIR, SLURM_GLOB))
fitsimage = CONTSUB_CUBE_FILE.format(imagebase=imagebase, sbid=sbid)

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

_linmos = cubestats_glob('cubeStats*linmos.contsub.txt')
if not _linmos:
    sys.exit("ERROR: no cubeStats*linmos.contsub.txt in " + str(cubestats_dirs()) + "\n       Run this script from inside the SBID directory.")
cubestat_linmos_contsub = _linmos[0] #mosaic contsub statistic

start_freq, end_freq = get_Frequency_Range(cubestat_linmos_contsub)
freq_range = str(start_freq)+'--'+str(end_freq)
delta_freq_range = round((end_freq - start_freq), 1)
theoretical_rms_mjy = (cal_Theoretical_RMS(float(n_ant), tobs, chan_width))*1000.
mad_rms, med_rms = cal_mosaic_Stats(cubestat_linmos_contsub)
n_bad_chan, mosaic_bad_chan, QC_badchan_id, QC_badchan_keyword = qc_Bad_Chans(cubestat_linmos_contsub, mad_rms, med_rms)
n_mdata_chan, mosaic_mdata_chan, QC_mdata_chan_id, QC_mdata_chan_keyword = qc_Missing_Data(cubestat_linmos_contsub)

# first_cat = get_FIRST(ra, dec)
# ---- Reference catalogue query - disabled for now; remove the "# [RACS] " prefixes to re-enable.
# [RACS] ref_cat_info = REF_CATALOGUES[REF_CATALOGUE]
# [RACS] if (QUERY_CATALOGUE or options.catalogue) and not options.no_catalogue:
# [RACS]     ref_cat = get_RefCatalogue(ra, dec)
# [RACS] else:
# [RACS]     ref_cat = ref_cat_info['outfile']
# [RACS]     if not os.path.isfile(fig_dir + '/' + ref_cat):
# [RACS]         print ("Catalogue query skipped", file=open(fig_dir + '/' + ref_cat, 'w'))

#Flagging statistic for spectral line
flagging_file = find_files(os.path.join(diagnostics_dir, FLAG_SUMMARY_GLOB))

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

#-------------------------------------------------------------------------------
# Diagnostic plots ordered as [contsub, image, residual]
#-------------------------------------------------------------------------------
def diag_by_kind(pattern):
    """
    Finding diagnostic plots and ordering them as [contsub, image, residual].
    A missing plot is replaced by a 'not available' placeholder, so the report still builds.
    """

    files = sorted(glob.glob(pattern)) if isinstance(pattern, str) else list(pattern)
    found = {}
    for f in files:
        name = os.path.basename(f).lower()
        kind = 'residual' if 'resid' in name else ('contsub' if 'contsub' in name else 'image')
        found.setdefault(kind, f)
    missing = [k for k in ('contsub', 'image') if k not in found]   # residual plots are not shown in the report
    if missing:
        print('WARNING: no %s plot matching %s - shown as "not available"' % (', '.join(missing),
              pattern if isinstance(pattern, str) else ('[' + ', '.join(os.path.basename(f) for f in files) + ']')))
    return [found.get(k) for k in ('contsub', 'image', 'residual')]


#-------------------------------------------------------------------------------
# Placeholder image for a missing plot
#-------------------------------------------------------------------------------
def make_NotAvailable(dir):
    """
    Placeholder image for missing diagnostic plots.
    """

    name = 'not_available.png'
    fig = plt.figure(figsize=(1.5, 1.5))
    fig.text(0.5, 0.5, 'not\navailable', ha='center', va='center', fontsize=14, color='gray')
    fig.savefig(dir + '/' + name, dpi=60)
    plt.close(fig)
    return dir + '/' + name


not_available_fig = make_NotAvailable(fig_dir)

cube_plots = diag_by_kind(cubestats_glob('*.cube*linmos*.png'))  #Mosaic statistic (all cubestats-* dirs)
#-------------------------------------------------------------------------------
# Find beamNoise / beamMinMax plots
#-------------------------------------------------------------------------------
def beam_diag_files(prefix):
    # beam-by-beam plots: diagnostics/ and every cubestats-* dir, any name containing the prefix
    # (residual plots can be named differently from the image/contsub ones, e.g. without '.cube')
    files = sorted(glob.glob(os.path.join(diagnostics_dir, '*' + prefix + '*.png'))) + cubestats_glob('*' + prefix + '*.png')
    print('%s plots found: %s' % (prefix, ', '.join(os.path.basename(f) for f in files) or 'none'))
    return files

beamNoise_plots = diag_by_kind(beam_diag_files('beamNoise')) #beam-by-beam statistic
beamMinMax_plots = diag_by_kind(beam_diag_files('beamMinMax')) #beam-by-beam statistic
spectra_glob = os.path.join(SOURCE_SPECTRA_DIR, SOURCE_SPECTRA_FILE.format(sbid=sbid, cid='*'))
source_spectra_fits = [os.path.join(SOURCE_SPECTRA_DIR, SOURCE_SPECTRA_FILE.format(sbid=sbid, cid=cid))
                       for cid in SOURCE_COMPONENTS]

thumb_cubeplots = []
thumb_beamNoise = []
thumb_beamMinMax = []
thumb_source_spectra = []

for plots, thumbs in [(cube_plots, thumb_cubeplots), (beamNoise_plots, thumb_beamNoise),
                      (beamMinMax_plots, thumb_beamMinMax)]:
    for k, image in enumerate(plots):
        if image is None:
            plots[k] = not_available_fig
            thumbs.append('not_available.png')
            continue
        thumb_img = 'thumb_'+ str(sizeX) + '_'+ os.path.basename(image)
        thumbs.append(thumb_img)
        make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

source_spectra_plots, source_spectra_labels = [], []
# Continuum image overlay
try:
    overlay_fig, overlay_plot = ContinuumOverlay_plot(sbid, field, tobs_hr)
except Exception as e:
    print('Could not make continuum overlay: ' + str(e))
    overlay_fig, overlay_plot = None, None
if overlay_fig:
    make_Thumbnail(overlay_fig, 'thumb_250_' + overlay_plot, 250, 250, fig_dir)

selavy_dict = read_Selavy(sbid)
for cid, specfile in zip(SOURCE_COMPONENTS, source_spectra_fits):
    if not os.path.isfile(specfile):
        print('Spectrum not found: ' + specfile)
        source_spectra_plots.append(None)
        source_spectra_labels.append('Component ' + cid + '<br>(not found)')
        thumb_source_spectra.append(None)
        continue
    try:
        noisefile = os.path.join(NOISE_SPECTRA_DIR, NOISE_SPECTRA_FILE.format(sbid=sbid, cid=cid))
        image, spec_plot, label = SourceSpectrum_plot(specfile, cid, selavy_dict, noisefile)
    except Exception as e:
        print('Could not plot spectrum ' + specfile + ': ' + str(e))
        source_spectra_plots.append(None)
        source_spectra_labels.append('Component ' + cid + '<br>(plot failed)')
        thumb_source_spectra.append(None)
        continue
    thumb_img = 'thumb_' + str(SPECTRA_THUMB_SIZE) + '_' + spec_plot
    make_Thumbnail(image, thumb_img, SPECTRA_THUMB_SIZE, SPECTRA_THUMB_SIZE, fig_dir)
    source_spectra_plots.append(image)
    source_spectra_labels.append(label)
    thumb_source_spectra.append(thumb_img)


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


# Noise spectra statistics (scripts 01-05)
noise_stats = None
noise_out = {}
try:
    os.makedirs(fig_dir + '/' + NOISE_DIR, exist_ok=True)
    ra_c, dec_c = get_FieldCentre(ra, dec)
    noise_stats = NoiseSpectra_stats(sbid, ra_c, dec_c, fig_dir + '/' + NOISE_DIR)
    if noise_stats is None:
        print('No noise spectra found - skipping noise statistics')
    else:
        selavy_tab = read_Selavy_table(sbid)
        fp_file = None
        for pattern in FOOTPRINT_GLOBS:
            fps = sorted(glob.glob(pattern.format(sbid=sbid)))
            if fps:
                fp_file = fps[0]
                break
        footprint = read_Footprint(fp_file) if fp_file else None
        noise_out['inout'] = NoiseInOut_plot(sbid, noise_stats)
        noise_out['dist'] = NoiseDist_plot(sbid, noise_stats)
        noise_out['sky'], noise_out['sky_summary'] = NoiseSky_plot(sbid, noise_stats, selavy_tab, ra_c, dec_c, footprint)
        noise_out['csky'], noise_out['csky_summary'] = NoiseSky_plot(sbid, noise_stats, selavy_tab, ra_c, dec_c,
                                                                     footprint, continuum=True)
        for key in ('inout', 'dist', 'sky', 'csky'):
            if noise_out.get(key):
                make_Thumbnail(fig_dir + '/' + noise_out[key], 'thumb_250_' + os.path.basename(noise_out[key]),
                               250, 250, fig_dir)
        for compno, pname in noise_stats['plots'].items():
            make_Thumbnail(fig_dir + '/' + pname, 'thumb_250_' + os.path.basename(pname), 250, 250, fig_dir)
except Exception as e:
    import traceback
    traceback.print_exc()
    print('Could not make noise statistics: ' + str(e))
    noise_stats = None


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
                    <title>%s</title>
                    """ % (css_style, REPORT_TITLE))

html.write("""<table width="100%" border="0">
                    <tr>
                        <th style="width:220"><img style="border:none;" src="" width="220" align="center"></th>
                        <th><h1 align="center"><font color="blue">""" + SURVEY_NAME + """</font> Spectral Line Data Validation Report</h1>
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
                                          dec.replace('.', ':', 2),   # mslist dd.mm.ss.ss -> dd:mm:ss.ss
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
                    <h2 align="middle">Beam Statistics</h2>
                    <table width="100%" border="1" style="table-layout:fixed">
                    <colgroup><col style="width:50.0000%"><col style="width:50.0000%"></colgroup>
                    <tr>
                    <th>Beam Image Cube</th>
                    <th>Continuum Subtracted Beam Cube</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1-percentile</p>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1-percentile</p>
                    </td>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{12}" target="_blank"><img src="{13}" width="{14}" height="{15}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM</p>
                    </td>
                    <td>
                    <a href="{16}" target="_blank"><img src="{17}" width="{18}" height="{19}" alt="thumbnail"></a>
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

# Catalog statistics: Selavy components and extracted source spectra
try:
    _cat_tab = read_Selavy_table(sbid)
except Exception as e:
    print('WARNING: could not read the Selavy catalogue for Catalog Statistics: %s' % e)
    _cat_tab = {}
_n_spec = len(glob.glob(os.path.join(SOURCE_SPECTRA_DIR, SOURCE_SPECTRA_FILE.format(sbid=sbid, cid='*'))))
if _cat_tab:
    _fint = np.array([v[4] for v in _cat_tab.values()])
    _n_bright = ('%d' % int(np.sum(_fint > CAT_BRIGHT_MJY))) if np.isfinite(_fint).any() else 'N/A (no flux_int column)'
    cat_stats_html = ('Number of components: %d<br>Number of components brighter than %g mJy: %s<br>'
                      'Total number of spectra extracted: %d'
                      % (len(_cat_tab), CAT_BRIGHT_MJY, _n_bright, _n_spec))
else:
    cat_stats_html = ('Number of components: N/A (no Selavy catalogue)<br>Number of components brighter than %g mJy: N/A<br>'
                      'Total number of spectra extracted: %d' % (CAT_BRIGHT_MJY, _n_spec))

html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Mosaic Statistics</h2>
                    <table width="100%" border="1" style="table-layout:fixed">
                    <colgroup><col style="width:16.6667%"><col style="width:16.6667%"><col style="width:16.6667%"><col style="width:16.6667%"><col style="width:16.6667%"><col style="width:16.6667%"></colgroup>
                    <tr>
                    <th colspan="3">Image Cube</th>
                    <th colspan="3">Continuum Subtracted Cube</th>
                    </tr>
		    <tr align="middle">
                    <td colspan="3">
                    <a href="{0}" target="_blank"><img src="{1}" width="{2}" height="{3}" alt="thumbnail"></a>
                    </td>
                    <td colspan="3">
                    <a href="{4}" target="_blank"><img src="{5}" width="{6}" height="{7}" alt="thumbnail"></a>
                    </td>
                    </tr>
                    <tr>
                    <th colspan="2">Number of Bad Quality Channels</th>
                    <th colspan="2">Number of Missing Data Channels</th>
                    <th colspan="2">Catalogue Statistics</th>
                    </tr>
                    <tr align="middle">
                    <td colspan="2"><br>{12}
                    <br><br><form action="{13}" method="get" target="_blank">
			<button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
                    </form>
		    </td>
                    <td colspan="2"><br>{14}
                    <br><br><form action="{15}" method="get" target="_blank">
			<button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
                    </form>
                    </td>
                    <td colspan="2" align="left">{16}
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
                               cat_stats_html))


html.write("""</td>
              </tr>
              </table>
              <h2 align="middle">Spectra Extracted from {0} Brightest Sources</h2>
              <table width="100%" border="1">
              """.format(N_SOURCE_SPECTRA))

if any(source_spectra_plots):
    for i, (image, thumb, label) in enumerate(zip(source_spectra_plots, thumb_source_spectra, source_spectra_labels)):
        if i % SPECTRA_NCOLS == 0:
            html.write('<tr align="middle">\n')
        if image is None:
            html.write('<td><br>' + label + '</td>\n')
        else:
            html.write("""<td>
              <a href="{0}" target="_blank"><img src="{1}" height="{2}" alt="thumbnail"></a>
              <br><p>{3}</p>
              </td>
              """.format(image, fig_dir + '/' + thumb, SPECTRA_THUMB_SIZE,
                         label))
        if i % SPECTRA_NCOLS == SPECTRA_NCOLS - 1 or i == len(source_spectra_plots) - 1:
            html.write('</tr>\n')
else:
    html.write('<tr align="middle"><td>No source spectra found (' + spectra_glob + ')</td></tr>\n')

html.write('</table>\n')

##### Noise spectra statistics
# Extra Validation Plots thumbnails are shown at the same width as the spectra thumbnails above
# (spectra are shown at height SPECTRA_THUMB_SIZE, so their width follows the figure aspect ratio)
EXTRA_THUMB_WIDTH = SPECTRA_THUMB_SIZE
for _t in thumb_source_spectra:
    if _t and os.path.isfile(fig_dir + '/' + _t):
        _w, _h = Image.open(fig_dir + '/' + _t).size
        EXTRA_THUMB_WIDTH = int(round(SPECTRA_THUMB_SIZE * _w / float(_h)))
        break

#-------------------------------------------------------------------------------
# HTML table cell with a thumbnail
#-------------------------------------------------------------------------------
def html_thumb(pname, caption, size=250, width=None):
    if not pname:
        return '<td>(not available)</td>\n'
    return ('<td><a href="{0}" target="_blank"><img src="{1}" width="{2}" alt="thumbnail"></a>'
            '<br>{3}</td>\n').format(fig_dir + '/' + pname,
                                             fig_dir + '/thumb_' + str(size) + '_' + os.path.basename(pname),
                                             width or EXTRA_THUMB_WIDTH, '<p>' + caption + '</p>' if caption else '')


html.write('<h2 align="middle">Spectral and Continuum Noise Statistics</h2>\n')
if noise_stats is None:
    html.write('<p>No noise spectra found (' + os.path.join(NOISE_SPECTRA_DIR, NOISE_SPECTRA_FILE.format(sbid=sbid, cid='*')) + ')</p>\n')
else:
    n3 = noise_stats['n3']
    comb = fig_dir + '/' + NOISE_DIR + '/SB%s_combined_all.txt' % sbid
    html.write(('<table width="100%" border="1" style="table-layout:fixed">\n<tr>'
               '<th>Spectral Noise: Inside vs Outside INNER_R deg<br>(Noise Spectra)</th>'
               '<th>RMS Noise vs Distance From Field Centre<br>(Noise Spectra)</th>'
               '<th>Spectral Noise Sky Distribution<br>(Noise Spectra)</th><th>Continuum Noise Sky Distribution<br>(Continuum Image)</th>'
               '<th>Distribution of Bright Components<br>(Continuum Image)</th></tr>\n'
               '<tr align="middle">\n').replace('INNER_R', '%g' % INNER_RADIUS_DEG))
    # No captions under the thumbnails (user request); only the combined-table button stays.
    inout_cell = html_thumb(noise_out.get('inout'), '')
    if os.path.isfile(comb):
        inout_cell = inout_cell.replace('</td>', '<form action="{0}" method="get" target="_blank">'
                                        '<button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>'
                                        '</form></td>'.format(comb))
    html.write(inout_cell)
    html.write(html_thumb(noise_out.get('dist'), ''))
    html.write(html_thumb(noise_out.get('sky'), ''))
    html.write(html_thumb(noise_out.get('csky'), ''))
    html.write(html_thumb(overlay_plot if overlay_fig else None, ''))
    html.write('</tr>\n</table>\n')

# ---- Reference catalogue (RACS-low) section - disabled for now; remove the "# [RACS] " prefixes to re-enable.
# [RACS] html.write("""
# [RACS]               <h2 align="middle">{2} sources within {1} &times; {1} sq degree ({3})</h2>
# [RACS]               <table>
# [RACS]               <tr align="middle">
# [RACS]               <td id="cell">
# [RACS]               <form action="{0}" method="get" target="_blank">
# [RACS]                  <button type = "submit" style="font-size:15px; width=50%; height=50%">Click here</button>
# [RACS]               </form>
# [RACS]               """.format(fig_dir+'/' + ref_cat, CAT_FIELD_WIDTH_DEG, ref_cat_info['label'], ref_cat_info['ref']))
# [RACS] html.write('</td>\n</tr>\n</table>\n')


##### Finish HTML report with generated time stamp
html.write("""


              <p align=left>
              Generated at {0} <br>
              <i> Report bugs to
              <a href="mailto:{1}">{2}</a></i>
              </p>
              </body>
              </html>""".format(str(datetime.now()), CONTACT_EMAIL, CONTACT_NAME))


html.close()
print ("Spectral line validation report written to '{0}'.".format(html_name))
