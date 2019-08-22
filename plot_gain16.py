# This script plots the gain solution, phase vs time and amplitude vs time for 12 antennas (XX and YY)
# It will read off the cont_gain*beam*.tab from ASKAPSoft that is stored in processed/M83_A.
# CASA browsetable doesn't convert imaginary and real part. Thus, this script is written. 
# Gain npol, nant, beam (always 0 as 1--3 array is impty), sampled time
# To run, casapy -c <script> <filename>

import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.pylab as pylab


infile=sys.argv[3]

#tb.open(infile)
tb.open(infile, nomodify=False)
tb.removerows([0])
times=tb.getcol('TIME')
gain=tb.getcol('GAIN')
tb.close()

beamno=str(infile[-6:-4])
sbname=str(infile[15:21])

# plot the phase vs time 

params={'xtick.labelsize': 'xx-small',
        'ytick.labelsize': 'xx-small'}

pylab.rcParams.update(params)

fig=plt.figure(figsize=(11.69,8.27))

for i in xrange(16):
    num=str(i)
    ax=plt.subplot(4,4,i+1)
    ax.tick_params(labelbottom='off')
    ax.plot(times,np.angle(gain[0,i,0,:],deg=True))
    ax.set_title('Antenna '+num+' XX', fontsize=8)
    yvalue=np.angle(gain[0,i,0,:],deg=True)
    ymin=np.min(yvalue)
    ymax=np.max(yvalue)
    ax.set_ylim(ymin-10,ymax+10)
    
fig.suptitle(sbname+' Beam'+beamno)    
fig.text(0.5, 0.04, 'Time [s]', ha='center', va='center')
fig.text(0.06, 0.5, 'Gain Phase [deg]', ha='center', va='center', rotation='vertical')
plt.savefig('gain_phase_ant.XX.beam'+beamno+'.pdf')
plt.close()

# plot the amplitude vs time (should be near 1, no amplitude calibration on last self-cal loop)

params={'xtick.labelsize': 'xx-small',
        'ytick.labelsize': 'xx-small'}

pylab.rcParams.update(params)

fig1=plt.figure(figsize=(11.69,8.27))
for y in xrange(16):
    num=str(y)
    ax1=plt.subplot(4,4,y+1)
    ax1.tick_params(labelbottom='off')
    ax1.plot(times,np.absolute(gain[0,y,0,:]))
    ax1.set_title('Antenna '+num+' XX', fontsize=8)
    ax1.set_ylim(0.5,1.5)

fig1.suptitle(sbname+' Beam'+beamno)    
fig1.text(0.5, 0.04, 'Time [s]', ha='center', va='center')
fig1.text(0.06, 0.5, 'Normalised Amplitude', ha='center', va='center', rotation='vertical')    
plt.savefig('gain_amplitude_ant.XX.beam'+beamno+'.pdf')
plt.close()

params={'xtick.labelsize': 'xx-small',
        'ytick.labelsize': 'xx-small'}

pylab.rcParams.update(params)

fig2=plt.figure(figsize=(11.69,8.27))
for z in xrange(16):
    num=str(z)
    ax2=plt.subplot(4,4,z+1)
    ax2.tick_params(labelbottom='off')
    ax2.plot(times,np.angle(gain[1,z,0,:],deg=True))
    ax2.set_title('Antenna '+num+' YY', fontsize=8)
    yvalue=np.angle(gain[1,z,0,:],deg=True)
    ymin=np.min(yvalue)
    ymax=np.max(yvalue)
    ax2.set_ylim(ymin-10,ymax+10)

fig2.suptitle(sbname+' Beam'+beamno)    
fig2.text(0.5, 0.04, 'Time [s]', ha='center', va='center')
fig2.text(0.06, 0.5, 'Gain Phase [deg]', ha='center', va='center', rotation='vertical')    
plt.savefig('gain_phase_ant.YY.beam'+beamno+'.pdf')
plt.close()


params={'xtick.labelsize': 'xx-small',
        'ytick.labelsize': 'xx-small'}

pylab.rcParams.update(params)

fig3=plt.figure(figsize=(11.69,8.27))
for zz in xrange(16):
    num=str(zz)
    ax3=plt.subplot(4,4,zz+1)
    ax3.tick_params(labelbottom='off')
    ax3.set_ylim(0.5,1.5)
    ax3.plot(times,np.absolute(gain[1,zz,0,:]))
    ax3.set_title('Antenna '+num+' YY', fontsize=8)

fig3.suptitle(sbname+' Beam'+beamno)    
fig3.text(0.5, 0.04, 'Time [s]', ha='center', va='center')
fig3.text(0.06, 0.5, 'Normalised Amplitude', ha='center', va='center', rotation='vertical')    
plt.savefig('gain_amplitude_ant.YY.beam'+beamno+'.pdf')
plt.close()
