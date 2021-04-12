################################################################################
# This script generates an ASKAP DINGO spectral line cube HTML report.
#
# Compatibility: Python version 3.x
#
# Required directories and files:
#
# -- metadata/mslist-scienceData*txt
# -- metadata/mslist-cal*txt
# -- metadata/mslist-*101.txt
# -- slurmOutput/*sh
# -- image.restored.i.SB<SBID>.cube.contsub.fits
# -- diagnostics/cubestats-<field>/*txt
# -- diagnostics/*png
# -- diagnostics/Flagging_Summaries/*SL.ms.flagSummary
# (below two are required only for continuum residaul tests)
# -- image.restored.i.SB_ID.cube.contsub.fits
# -- selavy-cont-image*restored/selavy-image*islands.xml'
#
# Output files: all saved in /validation_spectral/Figures directory.
#
# To run:
#
# python DINGO_specline_qa.py -s 10991 -c
#
# Originally Written for ASKAP Wallaby data validation by Bi-Qing For
# Version Modified on 24 March 2020
#
# Changed to work for ASKAP DINGO specral line data validation by Jonghwan Rhee
# Modified on 1 September 2020
# Changed to add continuum residual test part for ASKAP DINGO specral line data
# validation by Jonghwan Rhee (Modified on 25 September 2020)
# Updated on 12 April 2021 by Jonghwan Rhee
#
################################################################################

from astropy.stats import median_absolute_deviation
from astropy.utils.exceptions import AstropyWarning
from astropy.io.fits import getheader, getdata
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

from scipy.constants import k as k_B
from scipy.optimize import curve_fit
from scipy.stats import iqr
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.patches as mpatches
from matplotlib import cm

from datetime import datetime
from PIL import Image
from glob import glob
from math import sqrt, pi
import subprocess
import warnings
import sys
import os
import re
from argparse import ArgumentParser

################################################################################
# Functions for the main program
################################################################################


def footprint(name):
    """ Beam position offsets for ASKAP footprint

    Parameters
    ----------
    name: ASKAP footprint name (square_6x6, closepack36)

    Returns
    -------
    offsets: beam position offsets
    """
    if name == 'square_6x6':
        offsets = np.array([[-0.5, 0.5], [0.5, 0.5], [-0.5, -0.5], [0.5, -0.5], [-1.5, 1.5], [-0.5, 1.5],
                            [0.5, 1.5], [1.5, 1.5], [1.5, 0.5], [1.5, -0.5], [1.5, -1.5], [0.5, -1.5],
                            [-0.5, -1.5], [-1.5, -1.5], [-1.5, -0.5], [-1.5, 0.5], [-2.5, 2.5], [-1.5, 2.5],
                            [-0.5, 2.5], [0.5, 2.5], [1.5, 2.5], [2.5, 2.5], [2.5, 1.5], [2.5, 0.5],
                            [2.5, -0.5], [2.5, -1.5], [2.5, -2.5], [1.5, -2.5], [0.5, -2.5], [-0.5, -2.5],
                            [-1.5, -2.5], [-2.5, -2.5], [-2.5, -1.5], [-2.5, -0.5], [-2.5, 0.5], [-2.5, 1.5]]) * 1.0
        offsets[:, 0] *= -1.0

    elif name == 'closepack36':
        offsets = np.array([[-2.75, -2.16506], [-1.75, -2.16506], [-0.75, -2.16506], [0.25, -2.16506],
                            [1.25, -2.16506], [2.25, -2.16506], [-2.25, -1.29904], [-1.25, -1.29904],
                            [-0.25, -1.29904], [0.75, -1.29904], [1.75, -1.29904], [2.75, -1.29904],
                            [-2.75, -0.433013], [-1.75, -0.433013], [-0.75, -0.433013], [0.25, -0.433013],
                            [1.25, -0.433013], [2.25, -0.433013], [-2.25, 0.433013], [-1.25, 0.433013],
                            [-0.25, 0.433013], [0.75, 0.433013], [1.75, 0.433013], [2.75, 0.433013],
                            [-2.75, 1.29904], [-1.75, 1.29904], [-0.75, 1.29904], [0.25, 1.29904],
                            [1.25, 1.29904], [2.25, 1.29904], [-2.25, 2.16506], [-1.25, 2.16506],
                            [-0.25, 2.16506], [0.75, 2.16506], [1.75, 2.16506], [2.75, 2.16506]])
        offsets[:, 0] *= -1.0

    else:
        print('Your footprint is not available')
        sys.exit(1)

    return offsets


def gauss(x, *p):
    """ A Gaussian function.

    Parameters
    ----------
    x: data
    *p: list of Gussian parameters

    Returns
    -------
    f: Gaussian function
    """
    A, mu, sigma = p

    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))


def get_FitsHeader(fitsimage):
    """ Get basic information from the processed data
        For example, synthesized beam size
    """
    hdu = getheader(fitsimage)
    bmaj = round(float(hdu['BMAJ']) * 3600., 2)  # arcsec
    bmin = round(float(hdu['BMIN']) * 3600., 2)  # arcsec

    return bmaj, bmin


def get_beaminfo(beamlogs_file):
    """ Getting restored beam information in the middle of the channel

    Parameters
    ----------
    beamlogs_file: path to beamlogs file

    Returns
    -------
    bmaj, bmin: synthesized major and minor axis beam sizes
    """
    f = np.loadtxt(beamlogs_file)
    i_mid = int(f.shape[0] / 2 - 1)
    bmaj = f[i_mid][1]  # arcsec
    bmin = f[i_mid][2]  # arcsec

    return bmaj, bmin


def convert_dec_coords(casa_dec):
    """ Convert CASA Dec coordinate string to a proper format (delimiter '.' to ':')

    parameters
    ----------
    casa_dec: list of CASA Dec coordinates

    Returns
    -------
    coords_dec: list of dec coordinates with the delimiter of ':'
    """
    coords_dec = []
    for d in casa_dec:
        dec = d.split('.')
        dec[2:] = ['.'.join(dec[2:])]
        coords_dec.append(':'.join(dec))

    return coords_dec


def get_field_center(ra, dec):
    """ Calculate the mid position of all interleaves centers

    Parameters
    ----------
    ra: list of ra center position (string) for each interleave
    dec: list of dec ceter position (string) for each interleave

    Returns
    -------
    pos_cent: center position
    """
    ra_tmp = []
    dec_tmp = []
    for i in range(len(ra)):
        c = SkyCoord(ra[i], dec[i], unit=(u.hourangle, u.deg))
        ra_tmp.append(c.ra.value)
        dec_tmp.append(c.dec.value)

    pos_cent = [np.mean(ra_tmp), np.mean(dec_tmp)]

    return pos_cent


def get_HIPASS(pos_cen, dRA, dDec):
    """ Get the HIPASS sources (HICAT; Meyer et al. 2004)
        within a square coverage of the field through VizieR.

    Parameters
    ----------
    pos_cen: list of the field centre (RA, DEC)
    dRA: size of RA to search for HIPASS sources (in deg)
    dDec: size of Dec to search for HIPASS sources (in deg)

    Returns
    -------
    hipass_cat: name of the file containing HIPASS sources in the field coverage
    """
    from astroquery.vizier import Vizier

    print("Retrieving HIPASS sources from Vizier. Depending on server connection, this might take a while....")

    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=['HIPASS', '_RAJ2000', '_DEJ2000', 'RVsp', 'Speak', 'Sint', 'RMS', 'Qual'], catalog='VIII/73/hicat', timeout=1000)

    hipass_result = v.query_region(SkyCoord(ra=pos_cen[0], dec=pos_cen[1], unit=(u.deg, u.deg), frame='icrs'), width=[dRA * u.deg], height=[dDec * u.deg])

    hipass_cat = 'hipass.txt'
    #hipass_result['VIII/73/hicat'].write(fig_dir + '/' + hipass_cat, format='ascii.fixed_width', delimiter=' ')
    with open(fig_dir + '/' + hipass_cat, 'w') as f:
        print(hipass_result['VIII/73/hicat'], file=f)

    return hipass_cat


def get_version(param):
    """ Get the latest ASKAPsoft version used for the data reduction.

    Parameters
    ----------
    param: processing parameter file name

    Returns
    -------
    askapsoft: version of ASKAPsoft used
    """

    line = subprocess.check_output(['tail', '-5', param])  # Grab the last 5 lines
    str_line = line.decode('utf-8')
    askapsoft = re.findall('ASKAPsoft\ version\ [0-9].+', str_line)[0].split()[-1]

    return askapsoft


def get_Flagging(flagging_file):
    """ Get flagging statistics from flagSummary files

    Parameters
    ----------
    flagging_file: path (name) to the flagging surmmary files

    Returns
    -------
    data_flagged_pct: Flagged data fraction (%)
    """
    with open(flagging_file, 'r') as f:
        flag_info_list = f.readlines()

    for line in flag_info_list:
        if re.search('Total flagging', line):
            TOKS = line.split()
            total_flagged_pct = float(TOKS[-2])  # data+autocorrelation
            total_uv = float(TOKS[7])

    with open(flagging_file, 'r') as flag_infile:
        LINES = flag_infile.readlines()[:6]

    N_Rec = 'nRec'
    N_Chan = 'nChan'

    # Search for keywords in the file
    for i in range(len(LINES)):
        line = LINES[i]
        if line.find(N_Rec) >= 0:
            TOKS = line.split()
            n_Rec = float(TOKS[2])
        if line.find(N_Chan) >= 0:
            TOKS = line.split()
            n_Chan = float(TOKS[2])

    autocorr_flagged_pct = (36 * n_Rec * n_Chan / total_uv) * 100.0
    data_flagged_pct = round(total_flagged_pct - autocorr_flagged_pct, 2)

    return data_flagged_pct


def get_Metadata(metafile):
    """ Get basic information on observed field.
        obs date, total integration time, fields (name, ra, dec),
        # of antennas, total bandwidth

    Parameters
    ----------
    metafile: path to the metafile

    Returns
    -------
    info: dict of metadata containing obs information
    """
    # The number of 48-MHz correlator blocks (PILOT survey)
    nBlocks = 6

    # String patterns to match for extracting information
    duration = 'Total elapsed time = '
    obs_date = 'Observed from'
    fields = 'Fields: '
    code = 'Code'
    antennas = 'Antennas'
    totbw = 'TotBW'

    info = dict()
    with open(metafile, 'r') as f:
        mslist = f.readlines()
        for i, line in enumerate(mslist):
            # total integration time (in hr)
            if re.search(duration, line):
                info['tobs'] = round(float(line.split()[10]), 2)
                info['tobs_hr'] = round(float(line.split()[10])/3600., 2)

            # observation date and time
            if re.search(obs_date, line):
                info['start_obs_date'] = line.split()[6]
                info['end_obs_date'] = line.split()[8]

            # number of fields (2 for interleaving obs)
            if re.search(fields, line):
                info['nfields'] = int(line.split()[-1])

            # number of antennas
            if re.search(antennas, line):
                info['n_ant'] = int(line.split()[5])

            if re.search(code, line):
                for j in range(info['nfields']):
                    info['field_'+str(j)] = mslist[i+j+1].split()[4::]
                    info['field_'+str(j)].pop(-2)

            # total bandwidth in MHz
            if re.search(totbw, line):
                info['total_obs_bw'] = float(mslist[i+1].split()[10]) * nBlocks / 1000.
    return info


def get_Metadata_processed(metafile_science):
    """ Get basic information on processed data.
        number of channel, channel width, bandwidth and central frequency, band

    Parameters
    ----------
    metafile_science: path to the metafile

    Returns
    -------
    nchan: number of processed channels
    chan_width: channel width (in kHz)
    bw: bandwidth of processed data
    cfreq: central frequency of processed data
    band: ASKAP band (1/2)
    """
    with open(metafile_science, 'r') as f:
        for line in f.readlines():
            if re.search('Merged Window', line):
                nchan = line.split()[7]
                chan_width = float(line.split()[10])    # in kHz
                bw = float(line.split()[11]) / 1000.    # in MHz
                cfreq = line.split()[12]    # in MHz

    if float(cfreq) < 1150.:
        band = 1
    else:
        band = 2
    return nchan, chan_width, bw, cfreq, band


def get_frequency_range(cubestat_contsub):
    """ Get frequency range of the mosaic.

    Parameters
    ----------
    cubestat_contsub: path to the cube stats file of mosaic data

    Returns
    -------
    start_freq, end_freq: starting and ending frequency of processed data
    """

    line = subprocess.check_output(['sed', '-n', '3p', cubestat_contsub])
    TOKS = line.split()
    start_freq = round(float(TOKS[1]), 3)

    line = subprocess.check_output(['tail', '-1', cubestat_contsub])
    TOKS = line.split()
    end_freq = round(float(TOKS[1]), 3)

    return start_freq, end_freq


def make_Thumbnail(image, thumb_img, sizeX, sizeY, dir):
    """ Make a thumbnail image of input image with a size

    Parameters
    ----------
    image: name of input image
    thumb_img: name of thumbnail image
    sizeX: thumbnamil size in x axis
    sizeY: thumbnamil size in y axis
    dir: path to the output directory

    Returns
    -------
    thumb_img: name of thumbnail image
    """

    size = sizeX, sizeY
    im = Image.open(image)
    im.thumbnail(size)
    im.save(dir + '/' + thumb_img)

    return thumb_img


def cal_beam_MADMFD(infile):
    """ Calculate the Meian Absolute Deviation of max flux density of each beam.

    Parameters
    ----------
    infile: path to the cube statistic file

    Returns
    -------
    mad_maxfdensity: MAD of max flux density
    """
    data = np.loadtxt(infile)
    maxfdensity = data[:, 8]
    mad_maxfdensity = round(median_absolute_deviation(maxfdensity), 2)

    return mad_maxfdensity


def cal_beam_AvgRMS(infile):
    """ Calculate the average RMS of each beam.

    Parameters
    ----------
    infile: path to the cube statistic file

    Returns
    -------
    avg_rms: mean RMS
    """

    data = np.loadtxt(infile)
    rms = data[:, 3]
    avg_rms = round(np.nanmean(rms), 2)

    return avg_rms


def cal_beam_MedRMS(infile):
    """ Calculate the median RMS of each beam.

    Parameters
    ----------
    infile: path to the cube statistic file

    Returns
    -------
    med_rms: median RMS
    """

    data = np.loadtxt(infile)
    rms = data[:, 3]
    med_rms = round(np.nanmedian(rms), 2)

    return med_rms


def cal_mosaic_Stats(infile):
    """ Calculate MAD RMS and mean MADFM for the mosaic cube.

    Parameters
    ----------
    infile: path to the cube statistic file

    Returns
    -------
    mad_rms, med_madfm: MAD RMS, mean MADFM
    """
    data = np.loadtxt(infile)
    rms, madfm = data[:, 3], data[:, 5]

    med_madfm = np.nanmedian(madfm)
    mad_rms = round(median_absolute_deviation(rms), 2)

    return mad_rms, med_madfm


def qc_Bad_Chans(infile, med_madfm):
    """ Check for bad channels in the mosaic cube.

    Parameters
    ----------
    infile: path to the cube statistic file
    med_madfm: mean MADFM

    Returns
    -------
    mad_rms, med_madfm: MAD RMS, mean MADFM
    """
    BAD_CHAN = []
    BAD_FREQ = []

    with open(infile, 'r') as stat_file:
        LINES = stat_file.readlines()[2:]

    # Deviation from the med_madfm.
    # Need to check with larger sample of data to decide the best value.
    value = med_madfm + 0.4

    for i in range(len(LINES)):
        line = LINES[i]
        TOKS = line.split()
        chan = TOKS[0]
        freq = TOKS[1]
        madfm = float(TOKS[5])
        if madfm > value:
            BAD_CHAN.append(chan)
            BAD_FREQ.append(freq)

    if BAD_CHAN == []:
        BAD_CHAN.append('none')
        QC_badchan_id = 'good'
    else:
        QC_badchan_id = 'bad'

    mosaic_bad_chan = 'mosaic_badchans.txt'
    print('Channel  Frequency (MHz)', file=open(fig_dir + '/' + mosaic_bad_chan, 'w'))
    print('-------  ---------------', file=open(fig_dir + '/' + mosaic_bad_chan, 'a'))
    for c, f in zip(BAD_CHAN, BAD_FREQ):
        print(f'{c:>6} {f:>14}', file=open(fig_dir + '/' + mosaic_bad_chan, 'a'))

    n_bad_chan = len(BAD_CHAN)

    # Check if number of bad channel recorded is 1. If yes, check if is it a none keyword.
    # If yes, number of bad channel should be 0.

    if n_bad_chan == 1:
        with open(fig_dir + '/' + mosaic_bad_chan) as f:
            if 'none' in f.read():
                n_bad_chan = 0
                print('yes')

    return n_bad_chan, mosaic_bad_chan, QC_badchan_id


def cal_Theoretical_RMS(n_ant, tobs, chan_width, band):
    """
    Calculating the theoretical rms noise for ASKAP.
    Assuming natural weighting and not taking into account fraction of flagged data.
    For band 1, Tsys/eta = 75 K @ 1018.5 MHz is assumed based on ASKAP tsys measurements
    For band 2, Tsys/eta = 69 K @ 1367.5 MHz is assumed based on ASKAP tsys measurements
    """

    if band == 1:
        tsys = 52.6       # K
    elif band == 2:
        tsys = 48.3       # K
    else:
        print('Your band is not available')
        sys.exit(1)

    antdiam = 12    # m
    aper_eff = 0.7  # aperture efficiency
    coreff = 0.8    # correlator efficiency
    npol = 2.0      # Number of polarisation, npol = 2 for images in Stokes I, Q, U, or V

    anteff = pi * (antdiam / 2)**2. * aper_eff
    SEFD = 2. * k_B * 1e26 * tsys / anteff
    rms_jy = SEFD / (coreff * sqrt(npol * n_ant * (n_ant - 1) * chan_width * tobs))

    return rms_jy


def cal_Beam_ExpRMS(FLAGSTAT, t_rms):
    """ Calculate the theoretical RMS of individual beam
        by considering the flagged percentage
        Simple scaling with non-flagged data fraction

    Parameters
    ----------
    FLAGSTAT: file containing flagging data fraction for each beam
    t_rms: theoretical rms expected based on observation and ASKAP specs

    Returns
    -------
    BEAM_EXP_RMS: expected rms based on non-flagged data
    """
    BEAM_EXP_RMS = []

    for stat in FLAGSTAT:
        # 1/sqrt(non-flagged fraction) * theoretical rms in mJy
        beam_Exp_RMS = 1 / sqrt(float(1 - stat / 100)) * t_rms
        BEAM_EXP_RMS.append(beam_Exp_RMS)

    return BEAM_EXP_RMS


def qc_BeamLogs(field, band):
    """ Evaluate if the synthesized beam deviate from the tolerance range
        (+/-0.05 % of the nominal 30 arcsec which is currently set up).
        Note that bmaj and bmin for the first few channels are always zero
        due to barycentric correction.

    Parameters
    ----------
    field: name of an interleaving field
    """
    file_dir = diagnostics_dir + '/cubestats-' + field
    basename = '/beamlog.image.restored.' + imagebase + field
    if band == 1:
        tolerance = [20, 21]
    elif band == 2:
        tolerance = [30 - 30 * 0.006, 30 + 30 * 0.006]
    else:
        print('Your band is not available')
        sys.exit(1)

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


def FlagStat_plot(FLAGSTAT, field):
    """ Plotting and visualising flagging statistics of 36 beams.

    Parameters
    ----------
    FLAGSTAT: statistics of flagged data fraction
    field: name of the field (one interleaving)

    Returns
    -------
    saved_fig: path to the directroy of the figure
    plot_name: name of the figure
    """
    # set inputs (path to inputs and output file names)
    file_dir = diagnostics_dir + '/cubestats-' + field
    basename = '/cubeStats-image.restored.' + imagebase + field

    title = 'Flagged Data Fraction (' + field + ')'
    plot_name = 'FlagStat_' + field+ '.png'
    saved_fig = fig_dir + '/' + plot_name

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size': 10,
              'font.family': 'Times New Roman'}

    pylab.rcParams.update(params)

    # plot fraction of flagged data for each beam
    fig = plt.figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    offsets = footprint('square_6x6')

    for i in range(36):
        x0, y0 = offsets[i]
        plt.scatter([x0], [y0], s=1500, c=[FLAGSTAT[i]], cmap='tab20b', edgecolors='k', vmin=0, vmax=100)
        plt.text(x0, y0, '%d'%i, fontsize=12, va='center', ha='center')

    [i.set_linewidth(1.5) for i in ax.spines.values()]
    ax.tick_params(axis='both', which='both', bottom=False, top=False,
                    right=False, left=False, labelbottom=False, labelleft=False)
    ax.set_title(title)
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)

    cb = plt.colorbar(pad=0.02)
    cb.set_label('Percentage (%)')
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def Beam_ExpRMSplot(BEAM_EXPRMS, field):
    """ Plotting and visualising expected RMS of 36 beams.

    Parameters
    ----------
    BEAM_EXPRMS: expected RMS
    field: name of an interleaving field

    Returns
    -------
    saved_fig: path to the directroy of the figure
    plot_name: name of the figure
    """
    title = 'Expected RMS (' + field + ')'
    plot_name = 'Exp_RMS_' + field + '.png'
    saved_fig = fig_dir + '/' + plot_name

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size': 10,
              'font.family': 'Times New Roman'}

    pylab.rcParams.update(params)

    # plot fraction of flagged data for each beam
    fig = plt.figure(1)
    ax = fig.add_subplot(111, aspect='equal')

    offsets = footprint('square_6x6')
    VMIN = round(np.min(BEAM_EXPRMS), 1)
    VMAX = round(np.max(BEAM_EXPRMS), 1)

    for i in range(36):
        x0, y0 = offsets[i]

        plt.scatter([x0], [y0], s=1500, c=[BEAM_EXPRMS[i]], cmap='GnBu', edgecolors='black', vmin=VMIN, vmax=VMAX)
        plt.text(x0, y0, '%d'%i, fontsize=12, va='center', ha='center')

    [i.set_linewidth(1.5) for i in ax.spines.values()]
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    right=False, left=False, labelbottom=False, labelleft=False)
    ax.set_title(title)
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)

    cb = plt.colorbar(pad=0.02)
    cb.set_label('mJy / beam')
    plt.clim(VMIN, VMAX)
    cb.ax.tick_params(labelsize=10)

    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def BeamLogs_QCplot(list_beams_id_label, field):
    """ Plot Beam logs QC

    Parameters
    ----------
    list_beams_id_label: QC list of beams
    field: name of an interleaving field

    Returns
    -------
    saved_fig: path to the directroy of the figure
    plot_name: name of the figure
    """

    # set output file names
    plot_name = 'beamlogs_qc_SB' + sbid + '_' + field + '.png'
    saved_fig = fig_dir + '/' + plot_name

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size': 10,
              'font.family': 'Times New Roman'}

    pylab.rcParams.update(params)

    legend_dict = {'PASS': '#02b15d', 'FAIL': '#df2525'}
    patchList = []

    # plot QC signals for each beam
    fig = plt.figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    offsets = footprint('square_6x6')

    for i in range(36):
        x0, y0 = offsets[i]

        if (list_beams_id_label[i] == 'pass'):
            color_code = '#02b15d'
        else:
            color_code = '#df2525'

        plt.scatter(x0, y0, s=1500, c=color_code, edgecolors='k')
        plt.text(x0, y0, '%d'%i, fontsize=12, va='center', ha='center')

    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    [i.set_linewidth(1.5) for i in ax.spines.values()]
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)
    ax.tick_params(axis='both', which='both', bottom=False, top=False,
                    right=False, left=False, labelbottom=False, labelleft=False)
    ax.set_title('Beam Logs (' + field + ')')
    plt.legend(handles=patchList, ncol=2, loc='lower left', frameon=False, mode='expand', bbox_to_anchor=(0., 1.0, 1.0, 0.8), borderaxespad=0.)
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def BeamStat_plot(item, field):
    """ Plotting and visualising statistics of 36 beams.

    Parameters
    ----------
    item: Statistics to be visualised
    field: name of an interleaving field

    Returns
    -------
    saved_fig: path to the directroy of the figure
    plot_name: name of the figure

    """

    file_dir = diagnostics_dir + '/cubestats-' + field
    basename = '/cubeStats-image.restored.' + imagebase + field

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size': 10,
              'font.family': 'Times New Roman'}

    pylab.rcParams.update(params)

    if item == 'MADMFD':
        # This is a conservative cut off based on M83 field. 0.5 mJy/beam is more realistic.
        title = 'MAD Max Flux Density (' + field + ')'
        plot_name = 'beamStat_MADMFD_' + field + '.png'
        saved_fig = fig_dir + '/' + plot_name

    elif item == 'Avg_RMS':
        title = 'Mean RMS (' + field + ')'
        plot_name = 'beamStat_AvgRMS_' + field + '.png'
        saved_fig = fig_dir + '/' + plot_name

    elif item == 'Med_RMS':
        title = 'Median RMS (' + field + ')'
        plot_name = 'beamStat_MedRMS_' + field + '.png'
        saved_fig = fig_dir + '/' + plot_name

    # plot beam stats for each beam
    fig = plt.figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    offsets = footprint('square_6x6')

    beamstat_all = []
    for i in range(36):
        infile = file_dir + basename + '.beam%02d.txt' % (i)
        if os.path.isfile(infile):
            if item == 'MADMFD':
                beamstat = cal_beam_MADMFD(infile)
            elif item == 'Avg_RMS':
                beamstat = cal_beam_AvgRMS(infile)
            elif item == 'Med_RMS':
                beamstat = cal_beam_MedRMS(infile)

        else:
            print('Failed to read cubeStats field for beams')
            sys.exit(1)

        beamstat_all.append(beamstat)

    for i in range(36):
        x0, y0 = offsets[i]
        plt.scatter([x0], [y0], s=1500, c=[beamstat_all[i]],
                    cmap='GnBu', edgecolors='black', vmin=min(beamstat_all), vmax=max(beamstat_all))
        plt.text(x0, y0, '%d'%i, fontsize=12, va='center', ha='center')

    [i.set_linewidth(1.5) for i in ax.spines.values()]
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    right=False, left=False, labelbottom=False, labelleft=False)
    ax.set_title(title)
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)

    cb = plt.colorbar(pad=0.02)
    cb.ax.tick_params(labelsize=10)
    cb.set_label('mJy / beam')
    plt.clim(min(beamstat_all), max(beamstat_all))
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def qc_NoiseRank(spread):
    """ Evaluating the 1-percentile noise rank distribution.
    Sigma of the Gaussian function is used as a metric.
    """

    sigma = 0.5  # for a Gaussian

    if (spread <= 2.5 * sigma):
        qc_label = 'GOOD'
    elif (2.5 * sigma < spread < 4 * sigma):
        qc_label = 'OK'
    else:
        qc_label = 'BAD'

    return qc_label


def NoiseRank_histplot(nchan, field):
    """ Check distrubtion of 1 percentile noise

    Parameters
    ----------
    nchan: number of channels
    field: name of an interleaving field

    Returns
    -------
    saved_fig: path to the directroy of the figure
    plot_name: name of the figure
    ID_LABEL: list of 1-percentile noise quality labels for all beams

    """

    ID_LABEL = []
    plot_name = 'beam_1pctile_hist_SB' + sbid + '_' + field + '.png'
    saved_fig = fig_dir + '/' + plot_name
    file_dir = diagnostics_dir + '/cubestats-' + field
    basename = '/cubeStats-image.restored.' + imagebase + field

    params = {'axes.labelsize': 12,
              'axes.titlesize': 10,
              'xtick.labelsize': 9,
              'ytick.labelsize': 9,
              'font.size': 12,
              'font.family': 'Times New Roman'}

    pylab.rcParams.update(params)

    fig, axs = plt.subplots(6, 6, figsize=(16,9))
    fig.subplots_adjust(hspace=.004, wspace=.004, bottom=0.05)
    fig.text(0.5, -0.007, '1-percentile Noise Level (mJy / beam)', ha='center', va='center')
    fig.text(0.5, 1.02, '1-percentile Noise Distribution of ' + field, ha='center', va='center')
    fig.text(-0.001, 0.5, 'log N', ha='center', va='center', rotation='vertical')

    axs = axs.ravel()

    for i in range(36):
        infile = file_dir + basename + '.beam%02d.txt' % (i)
        data = np.loadtxt(infile)
        onepctile = data[:, 6]
        median_val = np.median(onepctile)
        # if statement is needed to rule out really bad data without having to do the Gaussian fitting
        if (median_val > 1000.0) or (median_val < -1000.0):
            ID_LABEL.append('bad')
            axs[i].set_xlim(-1, 1)
            axs[i].set_ylim(0, 3)
            axs[i].title.set_text('Beam%02d' % (i))
        else:
            upper_range = median_val + 3
            lower_range = median_val - 3
            x = onepctile[(onepctile < upper_range) & (onepctile > lower_range)]  # exclude outliers
            xmax_val = np.max(x)
            xmin_val = np.min(x)

            # Freedman-Diaconis rule. Nchan includes all processed channels, not excluding outliers.
            bin_width = 2 * iqr(x) * nchan**(-1 / 3)
            n_bins = int((xmax_val - xmin_val) / bin_width)

            hist, bins = np.histogram(
                onepctile, bins=n_bins, range=(xmin_val - 3, xmax_val + 3))
            with np.errstate(divide='ignore'):  # ignore division of zero
                N = np.log10(hist)   # get log N for y-axis
                N[N == -np.inf] = 0

            xcenter = (bins[:-1] + bins[1:]) / 2
            ymax_val = np.max(N)
            median_val_x = np.median(x)

            # Fitting a Gaussian and use spread (sigma) as a metric
            guess = [ymax_val, median_val_x, 5.0]
            coeff, var_matrix = curve_fit(gauss, xcenter, N, guess)
            spread = round(np.abs(coeff[2]), 1)
            ID_LABEL.append(qc_NoiseRank(spread))

            axs[i].bar(xcenter, N)
            axs[i].plot(xcenter, gauss(xcenter, *coeff), 'r-', lw=1)
            axs[i].set_xlim(xmin_val - 3, xmax_val + 3)
            axs[i].set_ylim(0, ymax_val + 3)
            axs[i].title.set_text('Beam%02d' % (i))

    plt.tight_layout()
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name, ID_LABEL


def NoiseRank_QCplot(list_id_label, field):

    params = {'axes.labelsize': 10,
              'axes.titlesize': 10,
              'font.size': 10,
              'font.family': 'Times New Roman'}

    pylab.rcParams.update(params)

    legend_dict = {'GOOD': '#02b15d', 'OK': '#f7c42a', 'BAD': '#df2525'}
    patchList = []

    plot_name = 'beam_1pctile_qc_SB' + sbid + '_' + field + '.png'
    saved_fig = fig_dir + '/' + plot_name


    # plot noise rank for each beam
    fig = plt.figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    offsets = footprint('square_6x6')

    for i in range(36):
        x0, y0 = offsets[i]
        if (list_id_label[i] == 'GOOD'):
            color_code = '#02b15d'
        elif (list_id_label[i] == 'OK'):
            color_code = '#f7c42a'
        else:
            color_code = '#df2525'

        plt.scatter([x0], [y0], s=1500, color=color_code, edgecolors='black')
        plt.text(x0, y0, '%d'%i, fontsize=12, va='center', ha='center')

    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    [i.set_linewidth(1.5) for i in ax.spines.values()]
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)

    plt.legend(handles=patchList, ncol=3, loc='lower left', frameon=False, mode='expand', bbox_to_anchor=(0., -0.08, 1.0, 0.8), borderaxespad=0.)
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    right=False, left=False, labelbottom=False, labelleft=False)
    plt.title('1-percentile Noise Rank (' + field + ')')
    plt.savefig(saved_fig, bbox_inches='tight')
    plt.close()

    return saved_fig, plot_name


def extract_spec(pos, cube):
    """ Extract a spectrum from a data cube

    parameters
    ----------
    pos: pixel positions of sources
    cube: data cube name where a spectrum is extracted

    Returns
    -------
    spec: array containing extracted spectra in mJy
    """

    spec = None
    for i in range(len(pos[0])):
        spec_tmp = getdata(cube)[:, 0, pos[1][i], pos[0][i]] * 1e3
        if spec is None:
            spec = spec_tmp
        else:
            spec = np.c_[spec, spec_tmp]

    return spec


def remove_outliers(x):
    """ Mask out outliers as nan value in the extracted spectra

    Parameters
    ----------
    x: an extracted spectra array (freq, source)

    Returns
    -------
    spec: the spectra array with outliers masked
    """
    spec = x.copy()
    lo_quartile = np.nanpercentile(spec, 10, axis=0)
    up_quartile = np.nanpercentile(spec, 90, axis=0)
    IQR = 1.5 * (up_quartile - lo_quartile)
    lo_bound = lo_quartile - IQR
    up_bound = up_quartile + IQR
    mask = (spec < lo_bound) | (spec > up_bound)
    spec[mask] = np.nan

    return spec


def select_sources(catalog, pos_cen, n_source=11):
    """ Select brightest continuum sources based on integrated flux within the
        central 5 x 5 sq. deg area.

    Parameters
    ----------
    catalog: continuum catalogue (sorted by source brightness)
    pos_cen: field centre coordinate
    n_source: number of top brightest source to select

    Returns
    -------
    catalog_new: new catalog based on the selection criteria
    """
    ra = catalog['ra_deg_cont']
    dec = catalog['dec_deg_cont']
    mask = (abs(ra.data - cent[0]) <= 2.5) & (abs(dec.data - cent[1]) <= 2.5)
    catalog_new = catalog[mask][:n_source]

    return catalog_new


def cal_dist(ra, dec, cent):
    """ Calculate the distance from the field center

    Parameters
    ----------
    ra: ra coordinates
    dec: dec coordinates
    cent: field center coordinats (ra, dec)

    Returns
    -------
    dist: distance from the field centre
    """
    c = SkyCoord(cent[0], cent[1], unit=(u.deg, u.deg))
    coords = SkyCoord(ra.data, dec.data, unit=(u.deg, u.deg))
    dist = c.separation(coords)

    return dist


def rebin_spec(x, y, bin_size, method):
    """ Rebinning a spectrum with a bin size

    Parameters
    ----------
    x : x axis (frequency)
    y : a spectrum
    bin_size : bin size in velocity or frequency
    method : re-binning method ('avg' : average / 'sum' : sum)

    Returns
    -------
    x_bins: rebinned x axis
    spec_rebinned : a rebinned spectrum with a bin size

    """
    x_bins = x[::bin_size]
    x_rebin = x[int(bin_size/2)::bin_size]
    mask_nan = np.isnan(y)
    y[mask_nan] = 0.0

    y_rebin = None
    for i in range(y.shape[1]):
        if y_rebin is None:
            if method == 'avg':
                y_rebin = np.histogram(x, bins=x_bins, weights=y[:,i])[0] / np.histogram(x[~mask_nan], bins=x_bins)[0]
            elif method == 'sum':
                y_rebin = np.histogram(x, bins=x_bins, weights=y[:,i])[0]
        else:
            if method == 'avg':
                y_tmp = np.histogram(x, bins=x_bins, weights=y[:,i])[0] / np.histogram(x[~mask_nan], bins=x_bins)[0]
                y_rebin = np.c_[y_rebin, y_tmp]
            elif method == 'sum':
                y_tmp = np.histogram(x, bins=x_bins, weights=y[:,i])[0]
                y_rebin = np.c_[y_rebin, y_tmp]
    if x_rebin.shape[0] > y_rebin.shape[0]:
        x_rebin = x_rebin[:-1]

    return x_rebin, y_rebin



################################################################################
# Main program where it calls all the functions
################################################################################

# ignore astropy warnings
warnings.simplefilter('ignore', AstropyWarning)

parser = ArgumentParser(description='Run DINGO validation and produce an HTML report')
parser.add_argument('-s','--sbid', dest='sbid', required='true', help='Science SBID', type=str)
parser.add_argument('-i','--imagebase', dest='imagebase', default='i.SB%s.cube', help='Base string for images [default=%default]', type=str)
parser.add_argument('-c','--contsubTest', dest='contsubTest', action="store_true", help="Whether to run the contsub test as well [default=%default]")
options = parser.parse_args()

# Switch to run contsub test
do_contsub_test = options.contsubTest

# Set file names and directories
diagnostics_dir = 'diagnostics'
fig_dir = 'Figures'
sbid = options.sbid
html_name = 'index.html'

imagebase=options.imagebase + '.'
imagebase=imagebase.replace('%s',sbid)

if not os.path.isdir(fig_dir):
    os.system('mkdir -p ' + fig_dir)

# Read required files
metafile = sorted(glob('metadata/mslist-*txt'))[0]
metafile_science = sorted(glob('metadata/mslist-scienceData*txt'))[0]
param_file = sorted(glob('slurmOutput/*.sh'))
beamlogs_file = sorted(glob('diagnostics/cubestats-*/beamlog*beam00.txt'))

# Check if there is more than one parameter input .sh file in the slurmOutput directory.
# If it does, select the latest one.
# If more than one version is used. Reporting the latest version of ASKAPSoft for final data reduction.

if len(param_file) > 1:
    param = param_file[-1]
else:
    param = param_file[0]

# Get observation information from metadata
info_metadata = get_Metadata(metafile)
askapsoft = get_version(param)
nchan, chan_width, bw, cfreq, band = get_Metadata_processed(metafile_science)
tobs_hr = info_metadata['tobs_hr']
chan_width_kHz = chan_width
field_names = [info_metadata['field_0'][1], info_metadata['field_1'][1]]

# Calculate on-source time for each interleaving field
dat_count = []
if info_metadata['nfields'] > 1:
    for i in range(info_metadata['nfields']):
        dat_count.append(float(info_metadata['field_' + str(i)][4]))

    dat_count = np.array(dat_count)
    t_int = dat_count / np.sum(dat_count) * tobs_hr
else:
    t_int.append(tobs_hr)

# mosaic contsub statistic
cubestat_linmos_contsub = sorted(glob(diagnostics_dir + '/cubestats-*/cubeStats*linmos.contsub.txt'))
cubestat_linmos_contsub_final = sorted(glob(diagnostics_dir + '/cubeStats*cube.contsub.txt'))

# get frequency information
start_freq, end_freq = get_frequency_range(cubestat_linmos_contsub[0])
delta_freq_range = end_freq - start_freq
bw_processed = end_freq - start_freq

# get beam, theoretical rms, observed rms, bad channels
beams = None
theoretical_rms = None

for i in range(info_metadata['nfields']):
    if beams is None:
        bmaj, bmin = get_beaminfo(beamlogs_file[i])
        beams = [(bmaj, bmin)]
        theoretical_rms = [(cal_Theoretical_RMS(info_metadata['n_ant'], t_int[i]*3600, chan_width*1e3, band)) * 1e3]
        ra_cen = [info_metadata['field_'+str(i)][2]]
        dec_cen = [info_metadata['field_'+str(i)][3]]

    else:
        beams.append(get_beaminfo(beamlogs_file[i]))
        theoretical_rms.append((cal_Theoretical_RMS(info_metadata['n_ant'], t_int[i]*3600, chan_width*1e3, band)) * 1e3)
        ra_cen.append(info_metadata['field_'+str(i)][2])
        dec_cen.append(info_metadata['field_'+str(i)][3])

# Check bad channels in the final mosaicked data cube
mad_rms, med_madfm = cal_mosaic_Stats(cubestat_linmos_contsub_final[0])
n_bad_chan, mosaic_bad_chan, QC_badchan_id = qc_Bad_Chans(cubestat_linmos_contsub_final[0], med_madfm)

# Calculate the center of all the interleaves
dec_cen = convert_dec_coords(dec_cen)
field_cen = get_field_center(ra_cen, dec_cen)
size_ra, size_dec = 5.5, 5.5
hipass_cat = get_HIPASS(field_cen, size_ra, size_dec)

# Flagging statistic for spectral line
flagging_file = sorted(glob('diagnostics/Flagging_Summaries/*SL.ms.flagSummary'))

FLAG_STAT_A = []
FLAG_STAT_B = []

for i, flag_file in enumerate(flagging_file):
    flag_stat = get_Flagging(flag_file)
    if i < 36:
        FLAG_STAT_A.append(flag_stat)
    else:
        FLAG_STAT_B.append(flag_stat)

BEAM_EXP_RMS_A = cal_Beam_ExpRMS(FLAG_STAT_A, theoretical_rms[0])
BEAM_EXP_RMS_B = cal_Beam_ExpRMS(FLAG_STAT_B, theoretical_rms[1])


################################################################################
# Run continuum residual test based on Sambit's code for this test
if do_contsub_test:

    # Input files
    selavy_file = glob('./selavy-cont-image*restored/selavy-image*islands.xml')[0]
    fitscube = glob('image.restored.' + imagebase + 'contsub.fits')[0]

    # Read selavy continuum catalogue
    selavy_cat = parse_single_table(selavy_file).to_table(use_names_over_ids=True)
    selavy_cat.sort('flux_int', reverse=True)

    # Select 11 brightest sources (flux_int > 200 mJ) within the central 5 x 5 sq. deg
    cent = field_cen
    selavy_bright = select_sources(selavy_cat, cent)

    ra = selavy_bright['ra_deg_cont']
    dec = selavy_bright['dec_deg_cont']
    n_comp = selavy_bright['n_components']
    f_int = selavy_bright['flux_int']
    dist = cal_dist(ra, dec, cent)

    # Add central pix info
    f_int_all = np.append(f_int.data.data, np.zeros(1))
    dist_all = np.append(dist.value, np.zeros(1))
    n_comp_all = np.append(n_comp.data.data, np.zeros(1))
    mask_multi = n_comp_all > 1

    # Get header info from the mosaicked data cube
    hd = getheader(fitscube)
    wcs = WCS(hd).sub(['celestial'])

    # Convert sky coordinates to pixel coordinates
    pos = np.vstack((ra.data.data, dec.data.data)).T
    pix = np.transpose(wcs.wcs_world2pix(pos, 0))
    pix = np.rint(pix).astype('int64')
    pos_cent = np.array([cent])
    pix_cent = np.rint(wcs.wcs_world2pix(pos_cent, 0)).astype('int64')

    # Add central pixel to the pixel array
    pix_all = np.c_[pix, pix_cent.T]

    # Extract spectra
    spec = extract_spec(pix_all, fitscube)

    # Make frequency array
    chan = np.arange(hd['NAXIS4'])
    freq = (hd['CRVAL4'] + chan * hd['CDELT4']) * 1e-6    # in MHz

    # Cut out the edges of spectra that have some artefacts
    spec = spec[256:6689,:]
    freq = freq[256:6689]

    # Mark out outliers from the extracted spectra
    spec2 = remove_outliers(spec)

    # Calculate mean and std of spectra
    spec_mean = np.nanmean(spec2, axis=0)
    spec_std = np.nanstd(spec2, axis=0)

    # Rebin spectra with increasing a bin size (pow of 4)
    bins = np.power(4, np.arange(0,6))

    freq_rebin = {}
    spec_rebin = {}
    spec_rebin_std = None
    for bin_size in np.power(4, np.arange(1,6)):
        freq_rebin[f'bin_size:{bin_size}'] = rebin_spec(freq, spec2, bin_size, 'sum')[0]
        spec_rebin[f'bin_size:{bin_size}'] = rebin_spec(freq, spec2, bin_size, 'sum')[1]
        if spec_rebin_std is None:
            spec_rebin_std = np.nanstd(rebin_spec(freq, spec2, bin_size, 'sum')[1], axis=0)
        else:
            spec_rebin_std = np.c_[spec_rebin_std, np.nanstd(rebin_spec(freq, spec2, bin_size, 'sum')[1], axis=0)]

    # Concatenate stats of orginal spectra into rebinned stats
    spec_stats = np.c_[spec_std, spec_rebin_std]


    # Plot an individual spectrum extracted from brigtest continuum sources
    # and central pixel
    plt.close('all')
    fig, axs = plt.subplots(6, 2, figsize=(10,10), constrained_layout=True)
    cmap = cm.rainbow(np.linspace(0, 1, 6))
    for i, ax in enumerate(axs.flat):
        ax.plot(freq, spec[:,i], color=cmap[0])
        [i.set_linewidth(1.2) for i in ax.spines.values()]
        ax.set_ylim(-15, 15)
        if (i % 2 == 0) and (i < 10):
            ax.set_ylabel('Flux (mJy/beam)')
            ax.set_title(f'{f_int[i]:.1f} mJy source, {dist[i]:.1f} from the center')
        elif i == 10:
            ax.set_ylabel('Flux (mJy/beam)')
            ax.set_xlabel('Frequency (MHz)')
            ax.set_title(f'{f_int[i]:.1f} mJy source, {dist[i]:.1f} from the center')
        elif i == 11:
            ax.set_xlabel('Frequency (MHz)')
            ax.set_title('Central Pixel')
        else:
            ax.set_title(f'{f_int[i]:.1f} mJy source, {dist[i]:.1f} from the center')

    plt.savefig(fig_dir + '/plot_spectra.png', bbox_inches='tight')


    # Plot mean and std of spectra vs. Continuum flux
    fig2 = plt.figure(2, figsize=(6.4, 10.4))
    fig2.subplots_adjust(top=0.94, bottom=0.12, hspace=0.2)
    ax1 = fig2.add_subplot(311)
    ax1.plot(f_int_all[-1], spec_mean[-1], '*', color='C0', ms=12, label='Mean (central pix)')
    ax1.plot(f_int_all[:-1], spec_mean[:-1], 'o', color='C0', label='Mean')
    ax1.plot(f_int_all[-1], spec_std[-1], '*', color='C1', ms=12, label='Std (central pix)')
    ax1.plot(f_int_all[:-1], spec_std[:-1], 'o', color='C1', label='Std')
    ax1.plot(f_int_all[mask_multi], spec_std[mask_multi], 'o', ms=13, mfc='None', mec='C2', label='Multi Comp')
    [i.set_linewidth(1.2) for i in ax1.spines.values()]
    ax1.set_xlabel('Continuum Source Flux (mJy)')
    ax1.set_ylabel('Flux of Spectra (mJy)')
    ax1.legend(ncol=3, bbox_to_anchor=(1.02, 1.25), loc='upper right')

    # Plot mean, Std vs. Distance from the field centre
    ax2 = fig2.add_subplot(312)
    ax2.plot(dist_all[-1], spec_mean[-1], '*', color='C0', ms=12, label='Mean (central pix)')
    ax2.plot(dist_all[:-1], spec_mean[:-1], 'o', color='C0', label='Mean')
    ax2.plot(dist_all[-1], spec_std[-1], '*', color='C1', ms=12, label='Std (central pix)')
    ax2.plot(dist_all[:-1], spec_std[:-1], 'o', color='C1', label='Std')
    ax2.plot(dist_all[mask_multi], spec_std[mask_multi], 'o', ms=13, mfc='None', mec='C2')
    [i.set_linewidth(1.2) for i in ax2.spines.values()]
    ax2.set_xlabel('Distance from the Center (deg)')
    ax2.set_ylabel('Flux of Spectra (mJy)')

    # Plot RMS vs No. chans summed over
    ax3 = fig2.add_subplot(313)
    cmap2 = cm.rainbow_r(np.linspace(0, 1, spec_stats.shape[0]))
    for i in range(spec_stats.shape[0]):
        if i < (spec_stats.shape[0]-1):
            ax3.plot(bins, spec_stats[i], 'o-', color=cmap2[i], label=f'{f_int_all[i]:.1f} mJy')
        else:
            ax3.plot(bins, spec_stats[i], 'o-', color='k', label='Central pix')
            ax3.plot(bins, spec_stats[i][0]*np.sqrt(bins), '--', color='k')

    ax3.set_xlabel('No. of Channels summed over')
    ax3.set_ylabel('RMS (mJy/beam)')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    [i.set_linewidth(1.2) for i in ax3.spines.values()]
    ax3.legend(ncol=4, bbox_to_anchor=(0.5, -0.34), loc='center')
    plt.savefig(fig_dir + '/plot_spectra_stats.png', bbox_inches='tight')


    # Plot individual spectra extracted from brigtest continuum sources
    # and central pixel and overlaid spectra summed with different bin sizes

    fig3, axs = plt.subplots(6, 2, figsize=(10,12), constrained_layout=False)
    fig3.subplots_adjust(top=0.95, left=0.08, right=0.95, bottom=0.08, hspace=0.35, wspace=0.15)
    for i, ax in enumerate(axs.flat):
        ax.plot(freq, spec[:,i], color=cmap[0], label=r'N$_{\rm sum}$=1')
        ax.plot(freq_rebin['bin_size:4'], spec_rebin['bin_size:4'][:,i]/np.sqrt(4), color=cmap[1],     label=r'N$_{\rm sum}$=4')
        ax.plot(freq_rebin['bin_size:16'], spec_rebin['bin_size:16'][:,i]/np.sqrt(16), color=cmap[2],     label=r'N$_{\rm sum}$=16')
        ax.plot(freq_rebin['bin_size:64'], spec_rebin['bin_size:64'][:,i]/np.sqrt(64), color=cmap[3],     label=r'N$_{\rm sum}$=64')
        ax.plot(freq_rebin['bin_size:256'], spec_rebin['bin_size:256'][:,i]/np.sqrt(256), color=cmap[4],     label=r'N$_{\rm sum}$=256')
        ax.plot(freq_rebin['bin_size:1024'], spec_rebin['bin_size:1024'][:,i]/np.sqrt(1024), color=cmap[5],     label=r'N$_{\rm sum}$=1024')
        y_min = np.nanmin(spec_rebin['bin_size:1024'][:,i])/np.sqrt(1024)
        y_max = np.nanmax(spec_rebin['bin_size:1024'][:,i])/np.sqrt(1024)
        y_lim = max(abs(y_min), abs(y_max)) + 5.0
        ax.set_ylim(-y_lim, y_lim)

        [i.set_linewidth(1.2) for i in ax.spines.values()]
        if (i % 2 == 0) and (i < 10):
            ax.set_ylabel('Flux (mJy/beam)')
            ax.set_title(f'{f_int[i]:.1f} mJy source, {dist[i]:.1f} from the center')
        elif i == 10:
            ax.set_ylabel('Flux (mJy/beam)')
            ax.set_xlabel('Frequency (MHz)')
            ax.set_title(f'{f_int[i]:.1f} mJy source, {dist[i]:.1f} from the center')
        elif i == 11:
            ax.set_xlabel('Frequency (MHz)')
            ax.set_title('Central Pixel')
            leg = ax.legend(ncol=6, bbox_to_anchor=(-0.05, -0.45), loc='center')
        else:
            ax.set_title(f'{f_int[i]:.1f} mJy source, {dist[i]:.1f} from the center')
    fig3.suptitle('Total flux over block of N channels scaled by sqrt(N)', y=0.99, fontsize=14)
    plt.savefig(fig_dir + '/plot_spectra_summed.png', bbox_inches='tight')
    plt.close('all')



################################################################################
# HTML related tasks
################################################################################

# Making thumbnail images
sizeX = 70
sizeY = 70

cube_plots = sorted(glob(diagnostics_dir + '/cubestats-G*/*linmos*.png'))  # Mosaic statistic
cube_plots_final = sorted(glob(diagnostics_dir + '/cubePlot*.png'))  # Final Mosaic statistic
beamNoise_plots = sorted(glob(diagnostics_dir + '/beamNoise*.png'))  # beam-by-beam statistic
beamMinMax_plots = sorted(glob(diagnostics_dir + '/beamMinMax*.png'))  # beam-by-beam statistic

thumb_cubeplots = []
thumb_cubeplots_final = []
thumb_beamNoise = []
thumb_beamMinMax = []

for image in cube_plots:
    thumb_img = 'thumb_' + str(sizeX) + '_' + image.split('/')[-1]
    thumb_cubeplots.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

for image in cube_plots_final:
    thumb_img = 'thumb_' + str(sizeX) + '_' + image.split('/')[-1]
    thumb_cubeplots_final.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

for image in beamNoise_plots:
    thumb_img = 'thumb_' + str(sizeX) + '_' + image.split('/')[-1]
    thumb_beamNoise.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

for image in beamMinMax_plots:
    thumb_img = 'thumb_' + str(sizeX) + '_' + image.split('/')[-1]
    thumb_beamMinMax.append(thumb_img)
    make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)

if do_contsub_test:
    contsub_plots = sorted(glob(fig_dir+'/plot_spectra*png'))
    thumb_contsub_tests = []
    for image in contsub_plots:
        thumb_img = 'thumb_' + str(sizeX) + '_' + image.split('/')[-1]
        thumb_contsub_tests.append(thumb_img)
        make_Thumbnail(image, thumb_img, sizeX, sizeY, fig_dir)


# Measured MAD of Maximum Flux Density of each beam for different interleaves
beam_MADMFD_fig_A, MADMFD_plot_A = BeamStat_plot('MADMFD', field_names[0])
thumb_img_A = 'thumb_' + str(sizeX) + '_' + MADMFD_plot_A
make_Thumbnail(beam_MADMFD_fig_A, thumb_img_A, sizeX, sizeY, fig_dir)

beam_MADMFD_fig_B, MADMFD_plot_B = BeamStat_plot('MADMFD', field_names[1])
thumb_img_B = 'thumb_' + str(sizeX) + '_' + MADMFD_plot_B
make_Thumbnail(beam_MADMFD_fig_B, thumb_img_B, sizeX, sizeY, fig_dir)

# Measured median RMS of each beam for different interleaves
beam_Med_RMS_fig_A, MedRMS_plot_A = BeamStat_plot('Med_RMS', field_names[0])
thumb_img_A = 'thumb_'+ str(sizeX) + '_'+ MedRMS_plot_A
make_Thumbnail(beam_Med_RMS_fig_A, thumb_img_A, sizeX, sizeY, fig_dir)

beam_Med_RMS_fig_B, MedRMS_plot_B = BeamStat_plot('Med_RMS', field_names[1])
thumb_img_B = 'thumb_'+ str(sizeX) + '_'+ MedRMS_plot_B
make_Thumbnail(beam_Med_RMS_fig_B, thumb_img_B, sizeX, sizeY, fig_dir)

# Expected RMS of each beam for different interleaves
beamExpRMS_fig_A, beamExpRMS_plot_A = Beam_ExpRMSplot(BEAM_EXP_RMS_A, field_names[0])
thumb_img_A = 'thumb_' + str(sizeX) + '_' + beamExpRMS_plot_A
make_Thumbnail(beamExpRMS_fig_A, thumb_img_A, sizeX, sizeY, fig_dir)

beamExpRMS_fig_B, beamExpRMS_plot_B = Beam_ExpRMSplot(BEAM_EXP_RMS_B, field_names[1])
thumb_img_B = 'thumb_' + str(sizeX) + '_' + beamExpRMS_plot_B
make_Thumbnail(beamExpRMS_fig_B, thumb_img_B, sizeX, sizeY, fig_dir)

# Measured 1 percentile noise of each beam for different interleaves
beam_1pctile_histfig_A, onepctile_plot_A, list_id_label_A = NoiseRank_histplot(float(nchan), field_names[0])
thumb_img_A = 'thumb_' + str(sizeX) + '_' + onepctile_plot_A
make_Thumbnail(beam_1pctile_histfig_A, thumb_img_A, sizeX, sizeY, fig_dir)

beam_1pctile_histfig_B, onepctile_plot_B, list_id_label_B = NoiseRank_histplot(float(nchan), field_names[1])
thumb_img_B = 'thumb_' + str(sizeX) + '_' + onepctile_plot_B
make_Thumbnail(beam_1pctile_histfig_B, thumb_img_B, sizeX, sizeY, fig_dir)

beam_1pctile_QCfig_A, onepctile_QCplot_A = NoiseRank_QCplot(list_id_label_A, field_names[0])
thumb_img_A = 'thumb_' + str(sizeX) + '_' + onepctile_QCplot_A
make_Thumbnail(beam_1pctile_QCfig_A, thumb_img_A, sizeX, sizeY, fig_dir)

beam_1pctile_QCfig_B, onepctile_QCplot_B = NoiseRank_QCplot(list_id_label_B, field_names[1])
thumb_img_B = 'thumb_' + str(sizeX) + '_' + onepctile_QCplot_B
make_Thumbnail(beam_1pctile_QCfig_B, thumb_img_B, sizeX, sizeY, fig_dir)

# Check synthesized beam size of each beam for different interleaves
list_beams_id_label_A = qc_BeamLogs(field_names[0], band)
BeamLogs_QCfig_A, BeamLogs_QCplot_A = BeamLogs_QCplot(list_beams_id_label_A, field_names[0])
thumb_img_A = 'thumb_' + str(sizeX) + '_' + BeamLogs_QCplot_A
make_Thumbnail(BeamLogs_QCfig_A, thumb_img_A, sizeX, sizeY, fig_dir)

list_beams_id_label_B = qc_BeamLogs(field_names[1], band)
BeamLogs_QCfig_B, BeamLogs_QCplot_B = BeamLogs_QCplot(list_beams_id_label_B, field_names[1])
thumb_img_B = 'thumb_' + str(sizeX) + '_' + BeamLogs_QCplot_B
make_Thumbnail(BeamLogs_QCfig_B, thumb_img_B, sizeX, sizeY, fig_dir)

# Check flagged data fraction of each beam for different interleaves
Flagged_fig_A, Flagged_plot_A = FlagStat_plot(FLAG_STAT_A, field_names[0])
thumb_img_A = 'thumb_' + str(sizeX) + '_' + Flagged_plot_A
make_Thumbnail(Flagged_fig_A, thumb_img_A, sizeX, sizeY, fig_dir)

Flagged_fig_B, Flagged_plot_B = FlagStat_plot(FLAG_STAT_B, field_names[1])
thumb_img_B = 'thumb_' + str(sizeX) + '_' + Flagged_plot_B
make_Thumbnail(Flagged_fig_B, thumb_img_B, sizeX, sizeY, fig_dir)



################################################################################
# Creating html report
################################################################################

css_style = """<style>
                  body {
                        background-color: white;
                        font-family: san-serif;
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
                         border-spacing: 0px;
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
                    <title>ASKAP DINGO Spectral Line Data Validation Report</title>
                    <h1 align="middle">DINGO Spectral Line Data Validation Report</h1>
                    """ % (css_style))

# Observations
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
                        <td>{4:.2f}</td>
                        <td>{5}</br>{6}</td>
                        <td>{7}</br>{8}</td>
                        <td>{9}</br>{10}</td>
                        <td>{11}""".format(sbid,
                                          info_metadata['n_ant'],
                                          info_metadata['start_obs_date'],
                                          info_metadata['end_obs_date'],
                                          info_metadata['tobs_hr'],
                                          field_names[0],
                                          field_names[1],
                                          ra_cen[0],
                                          ra_cen[1],
                                          dec_cen[0],
                                          dec_cen[1],
                                          info_metadata['total_obs_bw']))

# Processed data
html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Processed Data</h2>
                    <table width="100%" border="1">
                    <tr>
                        <th>ASKAPsoft<br>version*</th>
                        <th>Frequency Range<br>(MHz)</th>
                        <th>No. of Channels</th>
                        <th>Central Frequency<br>(MHz)</th>
                        <th>Channel Width<br>(kHz)</th>
                        <th>Synthesised Beam<br>(arcsec &times; arcsec)</th>
                        <th>On-source Time</br>(hr)</th>
                    </tr>
                    <tr align="middle">
                        <td>{0}</td>
                        <td>{1} &#8210; {2}</br>({3})</td>
                        <td>{4}</td>
                        <td>{5}</td>
                        <td>{6}</td>
                        <td>{7:.2f}&nbsp;&times;&nbsp;{8:.2f}&nbsp;(A)</br>{9:.2f}&nbsp;&times;&nbsp;{10:.2f}&nbsp;(B)</td>
                        <td>{11:.2f}&nbsp;(A)</br>{12:.2f}&nbsp;(B)</td>
                        """.format(askapsoft,
                                   str(start_freq),
                                   str(end_freq),
                                   round(bw_processed),
                                   nchan,
                                   cfreq,
                                   chan_width_kHz,
                                   beams[0][0],
                                   beams[0][1],
                                   beams[1][0],
                                   beams[1][1],
                                   t_int[0],
                                   t_int[1]))

# beam logs, flagging, rms expectation
html.write("""</td>
                    </tr>
                    </table>
                    <table width="100%" border="1">
                    <tr>
                        <th>Beam Logs (A/B)</th>
                        <th>Flagged Data Fraction (A/B)</th>
                        <th>Expected RMS (A/B)</th>
                    </tr>
                    <tr align="middle">
                        <td>
                        <a href="{2}" target="_blank"><img src="{3}" width="{0}" height="{1}" alt="thumbnail"></a>
                        <a href="{4}" target="_blank"><img src="{5}" width="{0}" height="{1}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{6}" target="_blank"><img src="{7}" width="{0}" height="{1}" alt="thumbnail"></a>
                        <a href="{8}" target="_blank"><img src="{9}" width="{0}" height="{1}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{10}" target="_blank"><img src="{11}" width="{0}" height="{1}" alt="thumbnail"></a>
                        <a href="{12}" target="_blank"><img src="{13}" width="{0}" height="{1}" alt="thumbnail"></a>
                        """.format(sizeX,
                                   sizeY,
                                   BeamLogs_QCfig_A,
                                   fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + BeamLogs_QCplot_A,
                                   BeamLogs_QCfig_B,
                                   fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + BeamLogs_QCplot_B,
                                   Flagged_fig_A,
                                   fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + Flagged_plot_A,
                                   Flagged_fig_B,
                                   fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + Flagged_plot_B,
                                   beamExpRMS_fig_A,
                                   fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + beamExpRMS_plot_A,
                                   beamExpRMS_fig_B,
                                   fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + beamExpRMS_plot_B))

# Image statistics (beam-based)
html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Image Statistics (beam-based) </h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Beam Image Cubes</th>
                    <th>Continuum Subtracted Beam Cubes</th>
                    <th>Residual Beam Cubes</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{2}" target="_blank"><img src="{3}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{4}" target="_blank"><img src="{5}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile (A/B)</p>
                    </td>
                    <td>
                    <a href="{6}" target="_blank"><img src="{7}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{8}" target="_blank"><img src="{9}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile (A/B)</p>
                    </td>
                    <td>
                    <a href="{10}" target="_blank"><img src="{11}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{12}" target="_blank"><img src="{13}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Min, Max, 1 percentile (A/B)</p>
                    </td>
                    <tr align="middle">
                    <td>
                    <a href="{14}" target="_blank"><img src="{15}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{16}" target="_blank"><img src="{17}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM (A/B)</p>
                    </td>
                    <td>
                    <a href="{18}" target="_blank"><img src="{19}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{20}" target="_blank"><img src="{21}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM (A/B)</p>
                    </td>
                    <td>
                    <a href="{22}" target="_blank"><img src="{23}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{24}" target="_blank"><img src="{25}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Stdev, MADFM (A/B)</p>
                    """.format(sizeX,
                               sizeY,
                               beamMinMax_plots[1],
                               fig_dir.split('/')[-1] + '/' + thumb_beamMinMax[1],
                               beamMinMax_plots[3],
                               fig_dir.split('/')[-1] + '/' + thumb_beamMinMax[3],
                               beamMinMax_plots[0],
                               fig_dir.split('/')[-1] + '/' + thumb_beamMinMax[0],
                               beamMinMax_plots[2],
                               fig_dir.split('/')[-1] + '/' + thumb_beamMinMax[2],
                               beamMinMax_plots[4],
                               fig_dir.split('/')[-1] + '/' + thumb_beamMinMax[4],
                               beamMinMax_plots[5],
                               fig_dir.split('/')[-1] + '/' + thumb_beamMinMax[5],
                               beamNoise_plots[1],
                               fig_dir.split('/')[-1] + '/' + thumb_beamNoise[1],
                               beamNoise_plots[3],
                               fig_dir.split('/')[-1] + '/' + thumb_beamNoise[3],
                               beamNoise_plots[0],
                               fig_dir.split('/')[-1] + '/' + thumb_beamNoise[0],
                               beamNoise_plots[2],
                               fig_dir.split('/')[-1] + '/' + thumb_beamNoise[2],
                               beamNoise_plots[4],
                               fig_dir.split('/')[-1] + '/' + thumb_beamNoise[4],
                               beamNoise_plots[5],
                               fig_dir.split('/')[-1] + '/' + thumb_beamNoise[5]))

# Continuum subtracted beam cubes
html.write("""</td>
                    </tr>
                    </table>
                    <table width="100%" border="1">
                    <tr>
                    <th colspan="4">Continuum Subtracted Beam Cubes</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{2}" target="_blank"><img src="{3}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{4}" target="_blank"><img src="{5}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>MAD Max Flux Density (A/B)</p>
                    </td>
                    <td>
                    <a href="{6}" target="_blank"><img src="{7}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{8}" target="_blank"><img src="{9}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>Median RMS Noise (A/B)</p>
                    </td>
                    <td>
                    <a href="{10}" target="_blank"><img src="{11}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{12}" target="_blank"><img src="{13}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{14}" target="_blank"><img src="{15}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{16}" target="_blank"><img src="{17}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <br><p>1-percentile Noise Rank (A/B)</p>
                    """.format(sizeX,
                               sizeY,
                               beam_MADMFD_fig_A,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + MADMFD_plot_A,
                               beam_MADMFD_fig_B,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + MADMFD_plot_B,
                               beam_Med_RMS_fig_A,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + MedRMS_plot_A,
                               beam_Med_RMS_fig_B,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + MedRMS_plot_B,
                               beam_1pctile_histfig_A,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + onepctile_plot_A,
                               beam_1pctile_QCfig_A,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + onepctile_QCplot_A,
                               beam_1pctile_histfig_B,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + onepctile_plot_B,
                               beam_1pctile_QCfig_B,
                               fig_dir.split('/')[-1] + '/thumb_' + str(sizeX) + '_' + onepctile_QCplot_B))

# Mosaic statistics
html.write("""</td>
                    </tr>
                    </table>
                    <h2 align="middle">Mosaic Statistics</h2>
                    <table width="100%" border="1">
                    <tr>
                    <th>Image Cube (A/B)</th>
                    <th>Continuum Subtracted Cube (A/B)</th>
                    <th>Residual Cube (A/B)</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{2}" target="_blank"><img src="{3}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{4}" target="_blank"><img src="{5}" width="{0}" height="{1}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{6}" target="_blank"><img src="{7}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{8}" target="_blank"><img src="{9}" width="{0}" height="{1}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{10}" target="_blank"><img src="{11}" width="{0}" height="{1}" alt="thumbnail"></a>
                    <a href="{12}" target="_blank"><img src="{13}" width="{0}" height="{1}" alt="thumbnail"></a>

                    """.format(sizeX,
                               sizeY,
                               cube_plots[1],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots[1],
                               cube_plots[4],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots[4],
                               cube_plots[0],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots[0],
                               cube_plots[3],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots[3],
                               cube_plots[2],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots[2],
                               cube_plots[5],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots[5]))

# Mosaic statistics (A/B combined)
html.write("""</td>
                    </tr>
                    </table>
                    <table width="100%" border="1">
                    <tr>
                    <th>Combined Image Cube</th>
                    <th>Combined Continuum Subtracted Cube</th>
                    <th>Combined Residual Cube</th>
                    <th>Number of Bad Channels</th>
                    </tr>
                    <tr align="middle">
                    <td>
                    <a href="{2}" target="_blank"><img src="{3}" width="{0}" height="{1}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{4}" target="_blank"><img src="{5}" width="{0}" height="{1}" alt="thumbnail"></a>
                    </td>
                    <td>
                    <a href="{6}" target="_blank"><img src="{7}" width="{0}" height="{1}" alt="thumbnail"></a>
                    </td>
                    <td id='{8}'>{9}
                    <form action="{10}" method="get" target="_blank">
                     <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
                    </form>
                    """.format(sizeX,
                               sizeY,
                               cube_plots_final[1],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots_final[1],
                               cube_plots_final[0],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots_final[0],
                               cube_plots_final[2],
                               fig_dir.split('/')[-1] + '/' + thumb_cubeplots_final[2],
                               QC_badchan_id,
                               n_bad_chan,
                               fig_dir.split('/')[-1] + '/' + mosaic_bad_chan))

# Continuum Residaul test for Mosaic cube
if do_contsub_test:
    html.write("""</td>
                        </tr>
                        </table>
                        <h2 align="middle">Continuum Residual Test***</h2>
                        <table width="100%" border="1">
                        <tr>
                        <th>Spectra from Brightest Sources and Field Centre</th>
                        <th>Statistics of Extracted Spectra and Spectra summed</th>
                        <th>Spectra and their total flux over N channel summed</th>
                        </tr>
                        <tr align="middle">
                        <td>
                        <a href="{2}" target="_blank"><img src="{3}" width="{0}" height="{1}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{4}" target="_blank"><img src="{5}" width="{0}" height="{1}" alt="thumbnail"></a>
                        </td>
                        <td>
                        <a href="{6}" target="_blank"><img src="{7}" width="{0}" height="{1}" alt="thumbnail"></a>
                        </td>
                        """.format(sizeX,
                                   sizeY,
                                   contsub_plots[0],
                                   fig_dir.split('/')[-1] + '/' + thumb_contsub_tests[0],
                                   contsub_plots[1],
                                   fig_dir.split('/')[-1] + '/' + thumb_contsub_tests[1],
                                   contsub_plots[2],
                                   fig_dir.split('/')[-1] + '/' + thumb_contsub_tests[2]))

# HIPASS soruces within the field observed
html.write("""</td>
              </tr>
              </table>
              <h2 align="middle">HIPASS sources within the central {0} &times; {1} sq degree**</h2>
              <table>
              <tr align="middle">
              <td id="cell">
              <form action="{2}" method="get" target="_blank">
                 <button type = "submit" style="font-size:20px; width=50%; height=50%">Click here</button>
              </form>
              """.format(size_ra, size_dec, fig_dir.split('/')[-1] + '/' + hipass_cat))


# Finish HTML report with generated time stamp
html.write("""
              </td>
              </tr>
              </table>

              <p align=left>
              * If more than one version of ASKAPsoft is used for the whole reduction, the latest one is reported.<br>
              ** Does not take into account field rotation.<br>
              *** Using spectra extracted from brightest continuum sources and field center positions.<br>
              Generated at {0} <br>
              <i> Report bugs to
              <a href="mailto:jonghwan.rhee@uwa.edu.au">Jonghwan Rhee</a></i>
              </p>
              </body>
              </html>""".format(str(datetime.now())))

html.close()
print("Spectral line validation report written to '{0}'".format(html_name))
