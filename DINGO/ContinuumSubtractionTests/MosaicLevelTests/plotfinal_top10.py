import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import array
import sys
import scipy
from matplotlib.font_manager import FontProperties
import PIL
from PIL import Image

#####################################################################

f0=1.2955e3			# MHz, for now: user defined
df=1.851851851845e-02		# MHz for now: user defined
nch=7776			# for now: user defined

#####################################################################

posta1=numpy.loadtxt('pspecs_top10_'+sys.argv[1]+'_postcontsub.dat')
postaa=posta1[posta1[:,3000]!=0,:]

postca=numpy.loadtxt('pcspec_top10_'+sys.argv[1]+'_postcontsub.dat')

contpars=numpy.loadtxt('contpars_top10_'+sys.argv[1]+'.dat')

#####################################################################

fint=contpars[posta1[:,3000]!=0,0]
dist=contpars[posta1[:,3000]!=0,1]
rad=contpars[posta1[:,3000]!=0,2]
decd=contpars[posta1[:,3000]!=0,3]

fqsa=numpy.zeros(nch)
for i in range(nch):
    fqsa[i]=f0+i*df

#####################################################################

posta=postaa[:,256:6689]
postc=postca[256:6689]
fqs=fqsa[256:6689]

#####################################################################

cmap = cm.rainbow(numpy.linspace(0,1,6))

tt=numpy.zeros(6)
for i in range(6):
    tt[i]=i
ttt=array.array('i',(0 for i in range(tt.shape[0])))
for i in range(tt.shape[0]):
    ttt[i]=int(tt[i])

chn=posta.shape[1]

fig, (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11) = plt.subplots(nrows=11, ncols=1, figsize=(6.4,46.9))
fig.subplots_adjust(hspace=1)

ax1.set_ylabel('flux (mJy/bm)')
ax1.set_xlabel('frequency (MHz)')
ax1.set_title(label="{0:.2f}".format(fint[0])+' mJy source')
#ax1.set_xlim(1298,1440)
ax1.plot(fqs[:],1000*posta[0,:],color=cmap[ttt[0]],label='N=1')

ax2.set_ylabel('flux (mJy/bm)')
ax2.set_xlabel('frequency (MHz)')
ax2.set_title(label="{0:.2f}".format(fint[1])+' mJy source')
ax2.plot(fqs[:],1000*posta[1,:],color=cmap[ttt[0]],label='N=1')

ax3.set_ylabel('flux (mJy/bm)')
ax3.set_xlabel('frequency (MHz)')
ax3.set_title(label="{0:.2f}".format(fint[2])+' mJy source')
ax3.plot(fqs[:],1000*posta[2,:],color=cmap[ttt[0]],label='N=1')

ax4.set_ylabel('flux (mJy/bm)')
ax4.set_xlabel('frequency (MHz)')
ax4.set_title(label="{0:.2f}".format(fint[3])+' mJy source')
ax4.plot(fqs[:],1000*posta[3,:],color=cmap[ttt[0]],label='N=1')

ax5.set_ylabel('flux (mJy/bm)')
ax5.set_xlabel('frequency (MHz)')
ax5.set_title(label="{0:.2f}".format(fint[4])+' mJy source')
ax5.plot(fqs[:],1000*posta[4,:],color=cmap[ttt[0]],label='N=1')

ax6.set_ylabel('flux (mJy/bm)')
ax6.set_xlabel('frequency (MHz)')
ax6.set_title(label="{0:.2f}".format(fint[5])+' mJy source')
ax6.plot(fqs[:],1000*posta[5,:],color=cmap[ttt[0]],label='N=1')

ax7.set_ylabel('flux (mJy/bm)')
ax7.set_xlabel('frequency (MHz)')
ax7.set_title(label="{0:.2f}".format(fint[6])+' mJy source')
ax7.plot(fqs[:],1000*posta[6,:],color=cmap[ttt[0]],label='N=1')

ax8.set_ylabel('flux (mJy/bm)')
ax8.set_xlabel('frequency (MHz)')
ax8.set_title(label="{0:.2f}".format(fint[7])+' mJy source')
ax8.plot(fqs[:],1000*posta[7,:],color=cmap[ttt[0]],label='N=1')

ax9.set_ylabel('flux (mJy/bm)')
ax9.set_xlabel('frequency (MHz)')
ax9.set_title(label="{0:.2f}".format(fint[8])+' mJy source')
ax9.plot(fqs[:],1000*posta[8,:],color=cmap[ttt[0]],label='N=1')

ax10.set_ylabel('flux (mJy/bm)')
ax10.set_xlabel('frequency (MHz)')
ax10.set_title(label="{0:.2f}".format(fint[9])+' mJy source')
ax10.plot(fqs[:],1000*posta[9,:],color=cmap[ttt[0]],label='N=1')

ax11.set_ylabel('flux (mJy/bm)')
ax11.set_xlabel('frequency (MHz)')
ax11.set_title('central pixel')
ax11.plot(fqs[:],1000*postc[:],color=cmap[ttt[0]],label='N=1')

fig.suptitle('Spectra at the positions of continuum sources')

plt.savefig(sys.argv[1]+'_top10_spectra.png',dpi=300,format='png')
plt.show()

#####################################################################

postcrms=numpy.nanstd(1000*postc[:])

postcmn=numpy.nanmean(1000*postc[:])

postrms=numpy.zeros(fint.shape[0])

for i in range (fint.shape[0]):
    postrms[i]=numpy.nanstd(1000*posta[i,:])

postmn=numpy.zeros(fint.shape[0])

for i in range (fint.shape[0]):
    postmn[i]=numpy.nanmean(1000*posta[i,:])

#####################################################################

fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
fig.subplots_adjust(hspace=0.5)

ax1.set_xlabel('continuum source flux (mJy)')
ax1.set_ylabel('mean / rms (mJy/bm)')
ax1.semilogx()
ax1.set_title('mean = open symbols, rms = filled symbols')
ax1.plot(fint,postrms,'bo',label='_nolegend_')
ax1.plot(min(fint)/2,postcrms,'cP',label='_nolegend_')
ax1.plot(fint,postmn,'bo',fillstyle="none",label='_nolegend_')
ax1.plot(min(fint)/2,postcmn,'cP',fillstyle="none",label='_nolegend_')

ax2.set_xlabel('distance of source from mosaic centre (deg)')
ax2.set_ylabel('mean / rms (mJy/bm)')
#ax2.semilogx()
ax2.plot(dist,postrms,'bo',label='_nolegend_')
ax2.plot(min(dist)/2,postcrms,'cP',label='_nolegend_')
ax2.plot(dist,postmn,'bo',fillstyle="none",label='_nolegend_')
ax2.plot(min(dist)/2,postcmn,'cP',fillstyle="none",label='_nolegend_')

fig.suptitle(sys.argv[1]+' sources after contsub')

plt.savefig(sys.argv[1]+'_top10_overview.png',dpi=300,format='png')
plt.show()

#####################################################################

apostcrms=numpy.zeros(6)

apostrms=numpy.zeros((posta.shape[0],6))

chn=posta.shape[1]

apostcrms[0]=numpy.nanstd(1000*postc[:])

for i in range(posta.shape[0]):
    apostrms[i,0]=numpy.nanstd(1000*posta[i,:])

for ii in range(5):
    apostc=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        apostc[jj]=numpy.nansum(postc[(pow(4,ii+1)*jj):(pow(4,ii+1)*(jj+1))])
    apostcrms[ii+1]=numpy.nanstd(1000*apostc[:])
    for kk in range(posta.shape[0]):
        aposta=numpy.zeros(chn/pow(4,ii+1))
        for jj in range(chn/pow(4,ii+1)):
            aposta[jj]=numpy.nansum(posta[kk,(pow(4,ii+1)*jj):(pow(4,ii+1)*(jj+1))])
        apostrms[kk,ii+1]=numpy.nanstd(1000*aposta[:])

chav=numpy.zeros(6)
for i in range (6):
  chav[i]=pow(4,i)

tvpost=numpy.zeros(6)
tvpost[:]=apostcrms[0]*numpy.sqrt(chav[:])

cmap = cm.rainbow(numpy.linspace(0,1,posta.shape[0]))

tt=numpy.zeros(posta.shape[0])
for i in range(posta.shape[0]):
    tt[i]=posta.shape[0]-i-1
ttt=array.array('i',(0 for i in range(tt.shape[0])))
for i in range(tt.shape[0]):
    ttt[i]=int(tt[i])

#####################################################################

fig=plt.figure()
fig.suptitle('RMS of total flux summed over blocks of channels')
plt.xlabel('no. of channels summed over')
plt.ylabel('RMS (mJy/bm)')
plt.semilogx()
plt.semilogy()
for i in range(posta.shape[0]):
    plt.plot(chav[:],apostrms[i,:],'o-',color=cmap[ttt[i]],label="{0:.2f}".format(fint[i])+' mJy')
plt.plot(chav[:],apostcrms[:],'ko-',label='central pixel')
plt.plot(chav[:],tvpost[:],'--',color='darkgrey',label='_nolegend')
plt.legend()
plt.savefig(sys.argv[1]+'_top10_chav.png',dpi=300,format='png')
plt.show()

#####################################################################

cmap = cm.rainbow(numpy.linspace(0,1,6))

tt=numpy.zeros(6)
for i in range(6):
    tt[i]=i
ttt=array.array('i',(0 for i in range(tt.shape[0])))
for i in range(tt.shape[0]):
    ttt[i]=int(tt[i])

chn=posta.shape[1]

fig, (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11) = plt.subplots(nrows=11, ncols=1, figsize=(6.4,46.9))
fig.subplots_adjust(hspace=1)

ax1.set_ylabel('flux (mJy/bm)')
ax1.set_xlabel('frequency (MHz)')
ax1.set_title(label="{0:.2f}".format(fint[0])+' mJy source')
ax1.plot(fqs[:],1000*posta[0,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[0,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax1.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax2.set_ylabel('flux (mJy/bm)')
ax2.set_xlabel('frequency (MHz)')
ax2.set_title(label="{0:.2f}".format(fint[1])+' mJy source')
ax2.plot(fqs[:],1000*posta[1,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[1,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax2.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))
ax2.legend()

ax3.set_ylabel('flux (mJy/bm)')
ax3.set_xlabel('frequency (MHz)')
ax3.set_title(label="{0:.2f}".format(fint[2])+' mJy source')
ax3.plot(fqs[:],1000*posta[2,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[2,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax3.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax4.set_ylabel('flux (mJy/bm)')
ax4.set_xlabel('frequency (MHz)')
ax4.set_title(label="{0:.2f}".format(fint[3])+' mJy source')
ax4.plot(fqs[:],1000*posta[3,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[3,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax4.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax5.set_ylabel('flux (mJy/bm)')
ax5.set_xlabel('frequency (MHz)')
ax5.set_title(label="{0:.2f}".format(fint[4])+' mJy source')
ax5.plot(fqs[:],1000*posta[4,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[4,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax5.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax6.set_ylabel('flux (mJy/bm)')
ax6.set_xlabel('frequency (MHz)')
ax6.set_title(label="{0:.2f}".format(fint[5])+' mJy source')
ax6.plot(fqs[:],1000*posta[5,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[5,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax6.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax7.set_ylabel('flux (mJy/bm)')
ax7.set_xlabel('frequency (MHz)')
ax7.set_title(label="{0:.2f}".format(fint[6])+' mJy source')
ax7.plot(fqs[:],1000*posta[6,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[6,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax7.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax8.set_ylabel('flux (mJy/bm)')
ax8.set_xlabel('frequency (MHz)')
ax8.set_title(label="{0:.2f}".format(fint[7])+' mJy source')
ax8.plot(fqs[:],1000*posta[7,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[7,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax8.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax9.set_ylabel('flux (mJy/bm)')
ax9.set_xlabel('frequency (MHz)')
ax9.set_title(label="{0:.2f}".format(fint[8])+' mJy source')
ax9.plot(fqs[:],1000*posta[8,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[8,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax9.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax10.set_ylabel('flux (mJy/bm)')
ax10.set_xlabel('frequency (MHz)')
ax10.set_title(label="{0:.2f}".format(fint[9])+' mJy source')
ax10.plot(fqs[:],1000*posta[9,:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(posta[9,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax10.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

ax11.set_ylabel('flux (mJy/bm)')
ax11.set_xlabel('frequency (MHz)')
ax11.set_title('central pixel')
ax11.plot(fqs[:],1000*postc[:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(postc[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax11.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))

fig.suptitle('Total flux over block of N channels scaled by sqrt(N)')

plt.savefig(sys.argv[1]+'_top10_spectrachav.png',dpi=300,format='png')
plt.show()

#####################################################################
