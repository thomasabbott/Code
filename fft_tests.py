# -*- coding: utf-8 -*-
"""
FFT channel position simulations, as an experiment for me
Created on Fri Nov 16 08:28:17 2018

@author: tabbott
"""

import numpy as np
from matplotlib import pyplot as plt

def testtone(ch,testch=30):
    dt = 1
    testch = 30
    tmax=2**8
    ftone = 1.0 * ch / tmax
    S = np.sin(2*np.pi*ftone*np.arange(tmax))# + 0.01 * np.random.randn(tmax)
    F = 20*np.log10(np.abs(np.fft.fft(S))[testch])
    return F


channel = np.arange(27,33,0.01)
R = [testtone(ch) for ch in channel]
R=R-np.max(R)


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(14,9), gridspec_kw=dict(height_ratios=[3],width_ratios=[1]))
tax = ax;
#tax.plot(freq,F,label='fft',color='b',lw=0.2,ms=10,marker='.');  tax.set_xlabel('frequency, Hz');
tax.plot(channel,R,label='channel %.1f' % ch,color='b',lw=1,ms=0);
tax.plot(channel+1,R,label='channel+1' % ch,color='k',lw=1,ms=0,alpha=0.3);
tax.plot(channel-1,R,label='channel-1' % ch,color='k',lw=1,ms=0,alpha=0.3);

tax.set_xlabel('Channel number of tone');
tax.legend(loc=0); 
tax.set_ylabel('response, dB'); tax.grid(); tax.set_title('FFT experiment - the power that comes out of channel 30')
tax.set_ylim([-60,0]);  tax.set_xlim([channel[0],channel[-1]]);  
