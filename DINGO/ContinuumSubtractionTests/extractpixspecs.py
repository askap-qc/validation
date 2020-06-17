import numpy
from astropy.io import fits
import sys
import math
import time
from astropy.wcs import WCS

f=numpy.loadtxt('selavy-results_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_gist.txt')

hdr=fits.getheader('image.restored.i.'+sys.argv[1]+'.cube.'+sys.argv[2]+'.'+sys.argv[3]+'.contsub.fits') 
wcs = WCS(hdr)

dat=fits.getdata('image.restored.i.'+sys.argv[1]+'.cube.'+sys.argv[2]+'.'+sys.argv[3]+'.contsub.fits')

dat1=fits.getdata('image.restored.i.'+sys.argv[1]+'.cube.'+sys.argv[2]+'.'+sys.argv[3]+'.fits')

ra0=hdr['crval1']
dec0=hdr['crval2']
st0=hdr['crval3']
fr0=hdr['crval4']
x0=hdr['crpix1']			
y0=hdr['crpix2']			
szx=hdr['naxis1']
szy=hdr['naxis2']
bmaj=hdr['bmaj']
bmin=hdr['bmin']
bpa=hdr['bpa']
dx=hdr['cdelt1']
dy=hdr['cdelt2']

dxy=(abs(dx)+abs(dy))/2
rad=numpy.ceil(abs(bmaj/dxy))		# both are in degrees
amaj=bmaj/dxy				# both are in degrees
amin=bmin/dxy				# both are in degrees
thet=math.radians(bpa)

# Initializing array to hold spectra, first channel for F_int (mJy)
specs=numpy.zeros((f.shape[0],dat.shape[0]))
specs1=numpy.zeros((f.shape[0],dat.shape[0]))

# Parameters of continuum sources
contpars=numpy.zeros((f.shape[0],2))

# Spectra of the central beams
cspec=dat[:,0,y0-1,x0-1]
cspec1=dat1[:,0,y0-1,x0-1]

for i in range(f.shape[0]):
    contpars[i,0]=f[i,2]
    ra=f[i,0]
    dec=f[i,1]
    rdiff=ra-ra0
    ddiff=dec-dec0
    contpars[i,1]=numpy.sqrt(pow(rdiff,2)+pow(ddiff,2))
#    x=int(round(x0+rdiff*math.cos(math.radians(y0))/dx-1))			# accounting for numpy array format
#    y=int(round(y0+ddiff/dy-1))			# accounting for numpy array format
    x=int(round(wcs.wcs_world2pix(ra,dec,st0,fr0,0)[0]))
    y=int(round(wcs.wcs_world2pix(ra,dec,st0,fr0,0)[1]))
    if (x>=0 and y>=0 and x<szx and y<szy):	# only sources within area
        specs[i,:]=dat[:,0,y,x]		# Jy/bm
        specs1[i,:]=dat1[:,0,y,x]	# Jy/bm
        
print time.localtime()			
numpy.savetxt('pspecs_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_postcontsub.dat',specs)
numpy.savetxt('pcspec_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_postcontsub.dat',cspec)
numpy.savetxt('pspecs_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_precontsub.dat',specs1)
numpy.savetxt('pcspec_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_precontsub.dat',cspec1)
numpy.savetxt('contpars_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'.dat',contpars)
