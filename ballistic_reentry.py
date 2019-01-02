# -*- coding: utf-8 -*-
"""
Ballistic reentry simulation, for interest

@author: tabbott
"""

import numpy as np
from matplotlib import pyplot as plt

global tstore, hstore, vstore, astore, machstore, cdstore, densitystore


def simulatelaunch(hinitial=120000.0,vinitial=0.0,area=3.8,mass=3175,ts=0.1,cd=0,verbose=2):
    
    global tstore, hstore, vstore, astore, machstore, cdstore, densitystore
    # constants
    g = 9.8                     # m.s^2
#    mass = 3175.0               # kg   Soyuz descent module is about 7000 lbs
#    area = np.pi*2.2/4.0        # m**2 of soyuz 2.2 m diameter of the blunt end
#    hinitial = 120000.0         # metres drop height
#   cd = 0 is code for cd=Soyuz.  Otherwise use cd from mach number lookup
    dragcoeff = np.array([[0.0, 0.65, 0.95, 1.1, 1.78, 2.5,  6.0, 20.0], \
                 [1.0, 1.0,  1.20, 1.4, 1.5,  1.55, 1.6,  1.6]])  # soyuz capsule - graphs found in a book
    atmosphere = np.genfromtxt('atmosphere.csv',delimiter=',',skip_header=1)  # US Standard atmosphere
    
    ts = 0.1       # timestep in seconds
    h=hinitial
    v=vinitial
    t=0.0
    a=-g
    
    tstore = np.empty(1)
    hstore = np.empty(1)
    vstore = np.empty(1)
    astore = np.empty(1)
    machstore = np.empty(1)
    cdstore = np.empty(1)
    densitystore = np.empty(1)
    powerstore = np.empty(1)
    
    # the simulation
    # height = up.  velocity = up.  accel = up.  force = up
    while h>0:
        density = np.interp(h,atmosphere[:,0],atmosphere[:,1])
        speedsound = np.interp(h,atmosphere[:,0],atmosphere[:,2])
        mach = np.abs(v/speedsound)
        if cd==0:
            cd = np.interp(mach,dragcoeff[0,:],dragcoeff[1,:])
        drag = - np.sign(v) * area * cd * (v**2.0 * density)/2.0
        power = drag*np.abs(v)
        force = drag -  mass*g
        a = force / mass
        v = v + a*ts
        h = h + v*ts
        t = t + ts
    
        tstore = np.vstack((tstore,t))
        hstore = np.vstack((hstore,h))
        vstore = np.vstack((vstore,v))
        astore = np.vstack((astore,a))
        machstore = np.vstack((machstore,mach))
        cdstore = np.vstack((cdstore,cd))
        densitystore = np.vstack((densitystore,density))
        powerstore = np.vstack((powerstore,power))
        
        if verbose>3:
            if np.floor(t/ts)%200 == 0:
                print 't=%4.0f  h=%7.0f  v=%7.0f  a=%7.1f' % (t,h,v,a)
        

    tstore=tstore[1:]; hstore=hstore[1:]; vstore=vstore[1:]; astore=astore[1:];
    machstore=machstore[1:]; cdstore=cdstore[1:]; densitystore=densitystore[1:]; powerstore=powerstore[1:];

    if verbose>0:
        print 'Run from %.0f km finished, mass=%.0f kg' % (hinitial/1000.0,mass)
        print 'Amax = %.1f g at h = %.1f km.  Vmax=%.1f m/s.  Vterminal = %.1f m/s.  Hmax = %.1f km'  \
                % ((np.max(astore)+g)/g,hstore[astore.argmax()]/1000.0,np.min(vstore),vstore[-1],np.max(hstore)/1000.0)

    
    if verbose>1:
        # plot results
        fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(14,9), \
                            gridspec_kw=dict(height_ratios=[3,3,3],width_ratios=[1,1]))
        # against time
        tax = ax[0,0]
        tax.plot(tstore,hstore/1000.0,label='height',color='b');
        tax.legend(loc=0); tax.set_xlabel('time, s'); tax.set_ylabel('height, km'); tax.grid()
        tax.set_title('Graphs against time.')
        
        tax = ax[1,0]
        tax.plot(tstore,vstore,label='velocity',color='c')
        tax.legend(loc=0); tax.set_xlabel('time, s'); tax.set_ylabel('velocity, m/s'); tax.grid()
        
        tax = ax[2,0]
        tax.plot(tstore,(astore+g)/g,label='acceleration, g_apparent',color='r')
        tax.legend(loc=0); tax.set_xlabel('time, s'); tax.set_ylabel('accel, g'); tax.grid()
        
        # against height
        tax = ax[0,1]
        tax.plot(hstore/1000.0,machstore,label='mach',color='m')
        tax.plot(hstore/1000.0,cdstore,label='cd',color='k')
        tax.plot(hstore/1000.0,np.log10(densitystore)+3,label='density log g/m3',color='b')
        tax.plot(hstore/1000.0,np.log10(powerstore+1.0)-3,label='drag power, log kW',color='r')
        tax.legend(loc=0); tax.set_xlabel('height, km'); tax.set_ylabel('mach and cd'); tax.grid()
        tax.set_title('Graphs against height. Hstart = %.1f km' % (hinitial/1000.0))
        #tax.set_xlim([0,60])
        
        tax = ax[1,1]
        tax.plot(hstore/1000.0,vstore,label='velocity',color='c')
        tax.plot(hstore[::50]/1000.0,vstore[::50],label='5-second ticks',color='c',lw=0,ms=6,marker='+')
        tax.legend(loc=0); tax.set_xlabel('height, km'); tax.set_ylabel('velocity, m/s'); tax.grid()
        tax.plot(hstore[vstore.argmin()]/1000.0,(np.min(vstore)),lw=0,ms=5,marker='o',color='c')
        tax.text(hstore[vstore.argmin()]/1000.0,(np.min(vstore))+300,'vmax=%.1f m/s' % (np.min(vstore)),color='c')
        #tax.set_xlim([0,60])
        
        tax = ax[2,1]
        tax.plot(hstore/1000.0,(astore+g)/g,label='acceleration, g_apparent',color='r')
        tax.plot(hstore[::50]/1000.0,(astore[::50]+g)/g,label='5-second ticks',color='r',lw=0,ms=6,marker='+')
        tax.plot(hstore[astore.argmax()]/1000.0,(np.max(astore)+g)/g,lw=0,ms=5,marker='o',color='r')
        tax.text(hstore[astore.argmax()]/1000.0+7,(np.max(astore)+g)/g-1,'amax=%.1f m/s^2' % ((np.max(astore)+g)/g),color='r')
        tax.legend(loc=0); tax.set_xlabel('height, km'); tax.set_ylabel('accel, g'); tax.grid()
        #tax.set_xlim([0,60])


# try some simulations
simulatelaunch(hinitial=50000,verbose=3,vinitial=+940)
#simulatelaunch(vinitial=2000.0,hinitial=1000.0,cd=0.02,mass=100,area=3.1*0.1*0.1)

