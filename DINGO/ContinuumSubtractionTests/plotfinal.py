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

prea1=numpy.loadtxt('pspecs_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_precontsub.dat')
prea=prea1[prea1[:,3000]!=0,:]

posta1=numpy.loadtxt('pspecs_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_postcontsub.dat')
posta=posta1[posta1[:,3000]!=0,:]

prec=numpy.loadtxt('pcspec_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_precontsub.dat')

postc=numpy.loadtxt('pcspec_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_postcontsub.dat')

contpars=numpy.loadtxt('contpars_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'.dat')

#####################################################################

fint=contpars[prea1[:,3000]!=0,0]*1000
dist=contpars[prea1[:,3000]!=0,1]

fqs=numpy.zeros(nch)
for i in range(nch):
    fqs[i]=f0+i*df

#####################################################################

fontP = FontProperties()
fontP.set_size('small')

fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(nrows=5, ncols=1, figsize=(6.4,9.6))
fig.subplots_adjust(hspace=1)

ax1.set_ylabel('flux (mJy/bm)')
ax1.set_title(label="{0:.2f}".format(fint[0])+' mJy source')
ax1.plot(fqs[:],1000*prea[0,:],color='r',label='before imcontsub')
ax1.plot(fqs[:],1000*posta[0,:],color='b',label='_nolegend_')
ax1.legend(prop=fontP)

ax2.set_ylabel('flux (mJy/bm)')
ax2.set_title(label="{0:.2f}".format(fint[1])+' mJy source')
ax2.plot(fqs[:],1000*prea[1,:],color='r',label='_nolegend_')
ax2.plot(fqs[:],1000*posta[1,:],color='b',label='after imcontsub')
ax2.legend(prop=fontP)

ax3.set_ylabel('flux (mJy/bm)')
ax3.set_title(label="{0:.2f}".format(fint[2])+' mJy source')
ax3.plot(fqs[:],1000*prea[2,:],color='r',label='before imcontsub')
ax3.plot(fqs[:],1000*posta[2,:],color='b',label='_nolegend_')
#ax3.legend(prop=fontP)

ax4.set_ylabel('flux (mJy/bm)')
ax4.set_xlabel('frequency (MHz)')
ax4.set_title(label="{0:.2f}".format(fint[3])+' mJy source')
ax4.plot(fqs[:],1000*prea[3,:],color='r',label='_nolegend_')
ax4.plot(fqs[:],1000*posta[3,:],color='b',label='after imcontsub')
#ax4.legend(prop=fontP)

ax5.set_ylabel('flux (mJy/bm)')
ax5.set_xlabel('frequency (MHz)')
ax5.set_title(label="{0:.2f}".format(fint[4])+' mJy source')
ax5.plot(fqs[:],1000*prea[4,:],color='r',label='_nolegend_')
ax5.plot(fqs[:],1000*posta[4,:],color='b',label='after imcontsub')
#ax5.legend(prop=fontP)

fig.suptitle('Spectra at brightest continuum source positions')

plt.savefig('finalplots/'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_1.png',dpi=300,format='png')
#plt.show()

#####################################################################

precrms=numpy.nanstd(1000*prec[:])
postcrms=numpy.nanstd(1000*postc[:])

precmn=numpy.nanmean(1000*prec[:])
postcmn=numpy.nanmean(1000*postc[:])

prerms=numpy.zeros(fint.shape[0])
postrms=numpy.zeros(fint.shape[0])

for i in range (fint.shape[0]):
    prerms[i]=numpy.nanstd(1000*prea[i,:])
    postrms[i]=numpy.nanstd(1000*posta[i,:])

premn=numpy.zeros(fint.shape[0])
postmn=numpy.zeros(fint.shape[0])

for i in range (fint.shape[0]):
    premn[i]=numpy.nanmean(1000*prea[i,:])
    postmn[i]=numpy.nanmean(1000*posta[i,:])

numpy.savetxt('postrms_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'.dat',postrms)

#####################################################################

fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
fig.subplots_adjust(hspace=0.5)

ax1.set_xlabel('continuum source flux (mJy)')
ax1.set_ylabel('mean / rms (mJy/bm)')
ax1.semilogx()
ax1.set_title('mean = open symbols, rms = filled symbols')
ax1.plot(fint,prerms,'ro',label='before imcontsub')
ax1.plot(fint,postrms,'bo',label='after imcontsub')
ax1.plot(min(fint)/2,precrms,'mP',label='_nolegend_') #label='central beam before imcontsub')
ax1.plot(min(fint)/2,postcrms,'cP',label='_nolegend_')
ax1.plot(fint,premn,'ro',fillstyle="none",label='_nolegend_')
ax1.plot(fint,postmn,'bo',fillstyle="none",label='_nolegend_')
ax1.plot(min(fint)/2,precmn,'mP',fillstyle="none",label='_nolegend_')
ax1.plot(min(fint)/2,postcmn,'cP',fillstyle="none",label='_nolegend_')
ax1.legend(prop=fontP)#,bbox_to_anchor=(1.2,0.4))

ax2.set_xlabel('distance of source from beam centre (deg)')
ax2.set_ylabel('mean / rms (mJy/bm)')
ax2.semilogx()
ax2.plot(dist,prerms,'ro',label='_nolegend_')
ax2.plot(dist,postrms,'bo',label='_nolegend_')
ax2.plot(min(dist)/2,precrms,'mP',label='beam centre, before')#label='central beam before imcontsub')
ax2.plot(min(dist)/2,postcrms,'cP',label='beam centre, after')
ax2.plot(dist,premn,'ro',fillstyle="none",label='_nolegend_')
ax2.plot(dist,postmn,'bo',fillstyle="none",label='_nolegend_')
ax2.plot(min(dist)/2,precmn,'mP',fillstyle="none",label='_nolegend_')
ax2.plot(min(dist)/2,postcmn,'cP',fillstyle="none",label='_nolegend_')
ax2.legend(prop=fontP)

fig.suptitle(sys.argv[1]+'   '+sys.argv[2]+'   '+sys.argv[3])#+', mean/rms = open/filled')

#plt.tight_layout()
plt.savefig('finalplots/'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_0.png',dpi=300,format='png')
#plt.show()

#####################################################################

sdst=numpy.sort(dist)

hdist=numpy.zeros(5)

for i in range(5):
    for j in range(dist.shape[0]):
        if dist[j]==sdst[sdst.shape[0]-1-i]:
            hdist[i]=j

fontP = FontProperties()
fontP.set_size('small')

fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(nrows=5, ncols=1, figsize=(6.4,9.6))
fig.subplots_adjust(hspace=1)

ax1.set_ylabel('flux (mJy/bm)')
ax1.set_title(label="{0:.3f}".format(dist[int(hdist[0])])+' deg from centre')
ax1.plot(fqs[:],1000*prea[int(hdist[0]),:],color='r',label='before imcontsub')
ax1.plot(fqs[:],1000*posta[int(hdist[0]),:],color='b',label='_nolegend_')
ax1.legend(prop=fontP)

ax2.set_ylabel('flux (mJy/bm)')
ax2.set_title(label="{0:.3f}".format(dist[int(hdist[1])])+' deg from centre')
ax2.plot(fqs[:],1000*prea[int(hdist[1]),:],color='r',label='_nolegend_')
ax2.plot(fqs[:],1000*posta[int(hdist[1]),:],color='b',label='after imcontsub')
ax2.legend(prop=fontP)

ax3.set_ylabel('flux (mJy/bm)')
ax3.set_title(label="{0:.3f}".format(dist[int(hdist[2])])+' deg from centre')
ax3.plot(fqs[:],1000*prea[int(hdist[2]),:],color='r',label='before imcontsub')
ax3.plot(fqs[:],1000*posta[int(hdist[2]),:],color='b',label='_nolegend_')
#ax3.legend(prop=fontP)

ax4.set_ylabel('flux (mJy/bm)')
ax4.set_xlabel('frequency (MHz)')
ax4.set_title(label="{0:.3f}".format(dist[int(hdist[3])])+' deg from centre')
ax4.plot(fqs[:],1000*prea[int(hdist[3]),:],color='r',label='_nolegend_')
ax4.plot(fqs[:],1000*posta[int(hdist[3]),:],color='b',label='after imcontsub')
#ax4.legend(prop=fontP)

ax5.set_ylabel('flux (mJy/bm)')
ax5.set_xlabel('frequency (MHz)')
ax5.set_title(label="{0:.3f}".format(dist[int(hdist[4])])+' deg from centre')
ax5.plot(fqs[:],1000*prea[int(hdist[4]),:],color='r',label='_nolegend_')
ax5.plot(fqs[:],1000*posta[int(hdist[4]),:],color='b',label='after imcontsub')
#ax5.legend(prop=fontP)

fig.suptitle('Spectra at continuum source positions furthest from centre')

plt.savefig('finalplots/'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_2.png',dpi=300,format='png')
#plt.show()

#####################################################################

aprecrms=numpy.zeros(6)
apostcrms=numpy.zeros(6)

aprerms=numpy.zeros((10,6))
apostrms=numpy.zeros((10,6))

chn=prea.shape[1]

aprecrms[0]=numpy.nanstd(1000*prec[:])
apostcrms[0]=numpy.nanstd(1000*postc[:])

for i in range(10):
    aprerms[i,0]=numpy.nanstd(1000*prea[i,:])
    apostrms[i,0]=numpy.nanstd(1000*posta[i,:])

for ii in range(5):
    aprec=numpy.zeros(chn/pow(4,ii+1))
    apostc=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aprec[jj]=numpy.nansum(prec[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        apostc[jj]=numpy.nansum(postc[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    aprecrms[ii+1]=numpy.nanstd(1000*aprec[:])
    apostcrms[ii+1]=numpy.nanstd(1000*apostc[:])
    for kk in range(10):
        aprea=numpy.zeros(chn/pow(4,ii+1))
        aposta=numpy.zeros(chn/pow(4,ii+1))
        for jj in range(chn/pow(4,ii+1)):
            aprea[jj]=numpy.nansum(prea[kk,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
            aposta[jj]=numpy.nansum(posta[kk,(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        aprerms[kk,ii+1]=numpy.nanstd(1000*aprea[:])
        apostrms[kk,ii+1]=numpy.nanstd(1000*aposta[:])

chav=numpy.zeros(6)
for i in range (6):
  chav[i]=pow(4,i)

cmap = cm.rainbow(numpy.linspace(0,1,10))

tt=numpy.zeros(10)
for i in range(10):
    tt[i]=10-i-1
ttt=array.array('i',(0 for i in range(tt.shape[0])))
for i in range(tt.shape[0]):
    ttt[i]=int(tt[i])

#####################################################################

fig=plt.figure()
fig.suptitle('RMS of total flux summed over blocks of channels')
plt.title('Before imcontsub')
plt.xlabel('no. of channels summed over')
plt.ylabel('RMS (mJy/bm)')
plt.semilogx()
plt.semilogy()
for i in range(10):
    plt.plot(chav[:],aprerms[i,:],'o-',color=cmap[ttt[i]],label="{0:.2f}".format(fint[i])+' mJy')
plt.plot(chav[:],aprecrms[:],'ko-',label='central pixel')
plt.legend()
plt.savefig('finalplots/'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_3.png',dpi=300,format='png')
#plt.show()

#####################################################################

fig=plt.figure()
fig.suptitle('RMS of total flux summed over blocks of channels')
plt.title('After imcontsub')
plt.xlabel('no. of channels summed over')
plt.ylabel('RMS (mJy/bm)')
plt.semilogx()
plt.semilogy()
for i in range(10):
    plt.plot(chav[:],apostrms[i,:],'o-',color=cmap[ttt[i]],label="{0:.2f}".format(fint[i])+' mJy')
plt.plot(chav[:],apostcrms[:],'ko-',label='central pixel')
plt.legend()
plt.savefig('finalplots/'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_4.png',dpi=300,format='png')
#plt.show()

#####################################################################

cmap = cm.rainbow(numpy.linspace(0,1,6))

tt=numpy.zeros(6)
for i in range(6):
    tt[i]=i
ttt=array.array('i',(0 for i in range(tt.shape[0])))
for i in range(tt.shape[0]):
    ttt[i]=int(tt[i])

chn=prea.shape[1]

fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows=4, ncols=1, figsize=(6.4,9.6))
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
#ax1.legend()

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
#ax3.legend()

ax4.set_ylabel('flux (mJy/bm)')
ax4.set_xlabel('frequency (MHz)')
ax4.set_title('central pixel')
ax4.plot(fqs[:],1000*postc[:],color=cmap[ttt[0]],label='N=1')
for ii in range(5):
    aposta=numpy.zeros(chn/pow(4,ii+1))
    afqs=numpy.zeros(chn/pow(4,ii+1))
    for jj in range(chn/pow(4,ii+1)):
        aposta[jj]=1000*numpy.nansum(postc[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
        afqs[jj]=numpy.nanmean(fqs[(pow(4,ii+1)*jj):(pow(4,ii+1)*jj+pow(4,ii+1)-1)])
    ax4.plot(afqs[:],aposta[:]/pow(2,ii+1),color=cmap[ttt[ii+1]],label='N='+str(pow(4,ii+1)))
#ax4.legend()

fig.suptitle('Total flux over block of N channels scaled by sqrt(N)')

plt.savefig('finalplots/'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_5.png',dpi=300,format='png')
#plt.show()

#####################################################################
