# WALLABY Spectral Line Validation

This script (wallaby_hi_val.py) generates an ASKAP WALLABY spectral line cube HTML report. The report contains basic observational information of the data cube, statistics and metrics for validating the data. 

## Requirements

- python 3.x
- astropy
- scipy
- matplotlib

## Required files and directories

- metadata/mslist-scienceData*txt
- metadata/mslist-cal*txt
- metadata/mslist-*101.txt
- slurmOutput/*sh
- image.restored.i.SB\<SBID\>.cube.contsub.fits (optional see comment)
- diagnostics/cubestats-\<field\>/*txt
- diagnostics/*png
- diagnostics/Flagging_Summaries/*SL.ms.flagSummary
- SpectralCube_BeamLogs/*.txt

## Usage

To run type: python \<script name\> \<SBID\>

## Output files

The html report is saved in the directory where the script is ran.  

## Description

## Useful links
- Further information: https://jira.csiro.au/browse/ACES-368 
- Useful resources from SPARCS: http://spacs.pbworks.com/w/page/126067640/dataquality 
- GASKAP validation: https://github.com/jd-au/gaskap-validation 
- Continuum substraction discussion page: https://confluence.csiro.au/display/askapsst/Continuum+Subtraction 
