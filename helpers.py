#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri June 3 17:18:26 2022

@author: rlm

helpers scripts for all other things
"""
#%


import datetime
import numpy as np
import scipy.signal as sig

def namestr(obj):
    try:
        return([name for name in globals() if globals()[name] is obj][0])
    except IndexError:
        return([name for name in locals() if locals()[name] is obj][0])

def dwnsmp(listofvals,n=None):
    rng = np.random.default_rng()
    if n is not None:
        if n<1: #fractional downsampling
            n = int(len(listofvals)*n)
        if n < len(listofvals):
            return(rng.choice(listofvals,size=n,replace=0))
        else:
            return(listofvals)


def printnow(inputstr=''):
    print(inputstr+'check: '+str(datetime.datetime.now()),flush=True)


rescalep = lambda vals: (vals-np.amin(vals))/(np.amax(vals)-np.amin(vals))

def dosmooth(vals,dtrndwind=50):
    if type(dtrndwind) is int:
        nsamp = dtrndwind
    elif type(dtrndwind) is float:
        nsamp = int(len(vals)*dtrndwind)
    smthsrc = sig.convolve(vals-np.mean(vals), np.ones(nsamp)/nsamp, 'same')
    return(smthsrc+np.mean(vals))


def dosmoothover(arr,param,smoothwind=50):
    smooth = np.zeros_like(arr[param])
    for ii in np.arange(arr.shape[0]):
        smooth[ii,:] = dosmooth(arr[param][ii,:],smoothwind)
    return(smooth)

def quickfft(times,values,posonly=1):
    dt = np.mean(np.diff(times))  # time step
    fs = 1 / dt  # sampling frequency

    # compute FFT
    fft_vals = np.fft.fft(values - np.mean(values))  # remove DC component
    freqs = np.fft.fftfreq(len(times), d=dt)

    if posonly:
        positive_freqs = freqs[freqs > 0]
        positive_fft = np.abs(fft_vals[freqs > 0])
        return(positive_freqs, positive_fft)
    else:
        return(freqs, fft_vals)
