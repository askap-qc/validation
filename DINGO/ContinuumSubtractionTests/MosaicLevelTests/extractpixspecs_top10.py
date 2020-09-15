import numpy
from astropy.io import fits
import sys
import math
import time
from astropy.wcs import WCS
import datetime

##################################
# Field centre - User defined !!!#
##################################

cra=130.90416805555552
cdec=-0.8448040740740742

# Reading in file and header

f=numpy.loadtxt('selavy-results_'+sys.argv[1]+'_gist.txt')
#f10=f[f[:,5]>=100,:]                              # continuum sources >100 mJy

hdr=fits.getheader('/scratch/twestmeier/image.restored.i.'+sys.argv[1]+'.cube.contsub.fits') 
wcs = WCS(hdr)

# Reading in important parameters from header

ra0=hdr['crval1']
dec0=hdr['crval2']
st0=hdr['crval3']
fr0=hdr['crval4']
x0=int(hdr['crpix1'])			
y0=int(hdr['crpix2'])
szx=hdr['naxis1']
szy=hdr['naxis2']
nch=hdr['naxis4']
bmaj=hdr['bmaj']
bmin=hdr['bmin']
bpa=hdr['bpa']
dx=hdr['cdelt1']
dy=hdr['cdelt2']

#######################################################
# Choosing all sources within a box - User defined !!!#
#######################################################

dra=2.5         # in degrees
ddec=2.5        # in degrees

dras=abs(f[:,0]-cra)

fa=f[dras[:]<=2.5,:]

ddecs=abs(fa[:,1]-cdec)

fb=fa[ddecs[:]<=2.5,:]

# Choosing top 10 sources

fsrti=numpy.sort(fb[:,5])
fsrt=fsrti[::-1]

f10=numpy.zeros((10,6))

for i in range(10):
    f10[i,:]=fb[fb[:,5]==fsrt[i],:]

# Initializing array to hold spectra

specs=numpy.zeros((f10.shape[0],nch))

# Parameters of continuum sources to be kept here

contpars=numpy.zeros((f10.shape[0],4))

# Extracting spectra at mosaic centre

x=int(round(float(wcs.wcs_world2pix(cra,cdec,st0,fr0,0)[0])))
y=int(round(float(wcs.wcs_world2pix(cra,cdec,st0,fr0,0)[1])))
cspec=fits.getdata('/scratch/twestmeier/image.restored.i.'+sys.argv[1]+'.cube.contsub.fits')[:,0,y,x]     # Jy/bm

# Extracting spectra at continuum source positions

for i in range(f10.shape[0]):
    contpars[i,0]=f10[i,5]
    ra=f10[i,0]
    dec=f10[i,1]
    rdiff=ra-cra
    ddiff=dec-cdec
    contpars[i,1]=numpy.sqrt(pow(rdiff,2)+pow(ddiff,2))       # distance from beam centre in degrees
    contpars[i,2]=rdiff
    contpars[i,3]=ddiff
    x=int(round(float(wcs.wcs_world2pix(ra,dec,st0,fr0,0)[0])))
    y=int(round(float(wcs.wcs_world2pix(ra,dec,st0,fr0,0)[1])))
    if (x>=0 and y>=0 and x<szx and y<szy):	# only sources within mosaic
        specs[i,:]=fits.getdata('/scratch/twestmeier/image.restored.i.'+sys.argv[1]+'.cube.contsub.fits')[:,0,y,x]		# Jy/bm

# Writing out files

numpy.savetxt('pcspec_top10_'+sys.argv[1]+'_postcontsub.dat',cspec)
numpy.savetxt('pspecs_top10_'+sys.argv[1]+'_postcontsub.dat',specs)
numpy.savetxt('contpars_top10_'+sys.argv[1]+'.dat',contpars)
