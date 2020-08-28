import numpy
from astropy.io import fits
import sys
import math
import time
from astropy.wcs import WCS

# Reading in files

f=numpy.loadtxt('selavy-results_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_gist.txt')

hdr=fits.getheader('image.restored.i.'+sys.argv[1]+'.cube.'+sys.argv[2]+'.'+sys.argv[3]+'.contsub.fits') 
wcs = WCS(hdr)

dat=fits.getdata('image.restored.i.'+sys.argv[1]+'.cube.'+sys.argv[2]+'.'+sys.argv[3]+'.contsub.fits')

dat1=fits.getdata('image.restored.i.'+sys.argv[1]+'.cube.'+sys.argv[2]+'.'+sys.argv[3]+'.fits')

# Reading in important parameters from header

ra0=hdr['crval1']
dec0=hdr['crval2']
st0=hdr['crval3']
fr0=hdr['crval4']
x0=hdr['crpix1']			
y0=hdr['crpix2']			
szx=hdr['naxis1']
szy=hdr['naxis2']

# Initializing array to hold spectra

specs=numpy.zeros((f.shape[0],dat.shape[0]))
specs1=numpy.zeros((f.shape[0],dat.shape[0]))

# Parameters of continuum sources to be kept here

contpars=numpy.zeros((f.shape[0],2))

# Extracting spectra at the bore sight

cspec=dat[:,0,y0-1,x0-1]
cspec1=dat1[:,0,y0-1,x0-1]

# Extracting spectra at continuum source positions

for i in range(f.shape[0]):
    contpars[i,0]=f[i,2]                                    # F_int (Jy)
    ra=f[i,0]
    dec=f[i,1]
    rdiff=ra-ra0
    ddiff=dec-dec0
    contpars[i,1]=numpy.sqrt(pow(rdiff,2)+pow(ddiff,2))     # distance from beam centre in degrees
    x=int(round(wcs.wcs_world2pix(ra,dec,st0,fr0,0)[0]))
    y=int(round(wcs.wcs_world2pix(ra,dec,st0,fr0,0)[1]))
    if (x>=0 and y>=0 and x<szx and y<szy):	# only sources within area covered by beam cube
        specs[i,:]=dat[:,0,y,x]		# Jy/bm
        specs1[i,:]=dat1[:,0,y,x]	# Jy/bm
        
# Writing out files	

numpy.savetxt('pspecs_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_postcontsub.dat',specs)
numpy.savetxt('pcspec_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_postcontsub.dat',cspec)
numpy.savetxt('pspecs_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_precontsub.dat',specs1)
numpy.savetxt('pcspec_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_precontsub.dat',cspec1)
numpy.savetxt('contpars_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'.dat',contpars)
