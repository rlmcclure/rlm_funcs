#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created January 2025

@author: rlmcclure
TODO
- generate quick spectra for 103 class
"""
#%% load packages 
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
# from astropy.io import fits
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from plotters import *
from helpers import *
#%%
def planck_law(wavelength, temp):
    """Calculate the blackbody flux using Planck's Law."""
    h = 6.626e-27  # Planck's constant in erg·s
    c = 3.0e10  # Speed of light in cm/s
    k = 1.38e-16  # Boltzmann constant in erg/K
    wavelength_cm = wavelength * 1e-8  # Convert Angstroms to cm
    return (2 * h * c**2) / (wavelength_cm**5 * (np.exp((h * c) / (wavelength_cm * k * temp)) - 1))

def generate_smooth_spectrum(wavelengths, temp=5000,normed=0):
    """Generate a smooth galaxy-like spectrum based on a blackbody curve."""
    peaklambda = 2.898e7 / temp
    if normed:
        return np.exp(-((wavelengths - peaklambda) / (peaklambda/2.8))**2) * (wavelengths / peaklambda)**-1.5
    else:
        return planck_law(wavelengths, temp) * (wavelengths / peaklambda)**-1.5

def add_absorption_feature(spectrum, wavelengths, center, width, depth):
    """Convolve the spectrum with an absorption line (Gaussian)."""
    absorption_profile = 1 - depth * np.exp(-0.5 * ((wavelengths - center) / width)**2)
    return spectrum * absorption_profile

def add_emission_feature(spectrum, wavelengths, center, width, relativestrength):
    """Add emission line as Gaussian profiles."""
    strength = relativestrength*np.amax(spectrum)
    emission_profile = strength * np.exp(-0.5 * ((wavelengths - center) / width)**2)
    spectrum+= emission_profile
    return spectrum
#%%
# Define wavelength range
wavelengths = np.linspace(3000, 9000, 1000)  # Angstroms

# Generate smooth continuum spectrum
# for wl in [3000,4500,5500,7000]:
normed=0;wabs=1;wemis=1;summed=1;cumul=1
abs = [(3933,20,.1),(3968,20,.1),(5896,4,.15),(5890,24,.153),(5270,20,.2),(5328,20,.2),(6800,80,.3),(7200,30,.22),(7772,1,.15),(7774,1,.15),(7775,1,.15),(8500,.5,.1)]

Tblue=11000;Tmid=8000;Tsunlike=6500;Tlow=3000
#blue stars

f = wbkg((10,5))#initialzie figure

spectrum_blue_cont = generate_smooth_spectrum(wavelengths,Tblue,normed)#+6
spectrum_blue = np.copy(spectrum_blue_cont)
if wabs:
    for newline,wid,dep in abs:
        spectrum_blue = add_absorption_feature(spectrum_blue, wavelengths, center=newline, width=wid, depth=dep)

    
plt.plot(wavelengths, spectrum_blue, label='Continuum, Blue, High Mass Stars T=%iK'%Tblue,c='cyan')#%3000)
#intermediate stars
spectrum_mid_cont = generate_smooth_spectrum(wavelengths,Tmid,normed)#+3
spectrum_mid = np.copy(spectrum_mid_cont)
if wabs:
    for newline,wid,dep in abs:
        newline *= 1.1
        spectrum_mid = add_absorption_feature(spectrum_mid, wavelengths, center=newline, width=wid, depth=dep)
plt.plot(wavelengths, spectrum_mid, label='Continuum, Intermediate Mass Stars T=%iK'%Tmid,c='yellow')#%4500)
#sun-like stars
spectrum_sunlike_cont = generate_smooth_spectrum(wavelengths,Tsunlike,normed)#+1
spectrum_sunlike = np.copy(spectrum_sunlike_cont)
if wabs:
    for newline,wid,dep in abs:
        newline *= .8
        spectrum_sunlike = add_absorption_feature(spectrum_sunlike, wavelengths, center=newline, width=wid, depth=dep)

plt.plot(wavelengths, spectrum_sunlike, label='Continuum, Sun-like Stars T=%iK'%Tsunlike,c='orange')#%5500)
#red stars
spectrum_low_cont = generate_smooth_spectrum(wavelengths,Tlow,normed)
spectrum_low = np.copy(spectrum_low_cont)
if wabs:
    for newline,wid,dep in abs:
        newline *= 1.1
        spectrum_low = add_absorption_feature(spectrum_low, wavelengths, center=newline, width=wid, depth=dep)

plt.plot(wavelengths, spectrum_low, label='Continuum, Low Mass Stars T=%iK'%Tlow,c='red')#%7000)
plt.legend(loc=1)
plt.grid(alpha=.8)
if normed:
    plt.ylim(0,1.7)
    plt.ylabel("Flux, normalized (Arbitrary Units)")
else:
    plt.ylabel("Flux, physically scaled")
if wabs:
    plt.title('with random absorption features')
plt.xlabel("Wavelength (Angstroms)")
# plt.yscale('log')
plt.show()

#%% build galaxy 
NLow = 10000
NSunlike = 1000
NInt = 100
NBlue = 10
f,ax = wbkg((10,5),rtax=1) #initalize figure
galsepect = NLow*spectrum_low+NSunlike*spectrum_sunlike+NInt*spectrum_mid+NBlue*spectrum_blue# Add emission lines (e.g., H-alpha at 6563, [OIII] at 5007 Angstroms)
emis = [(6563,3,.3),(5007,1,.1)]
if wemis:
    for newline,wid,dep in emis:
        galsepect = add_emission_feature(galsepect, wavelengths, center=newline, width=wid, relativestrength=dep)
ax2 = ax.twinx()
if summed:
    ax.plot(wavelengths, NBlue*spectrum_blue, label='N=%i, Blue, High Mass Stars T=%iK, total'%(NBlue,Tblue),c='cyan')#%3000)
    ax.plot(wavelengths, NBlue*spectrum_blue_cont, ls=':',c='cyan',lw=.5)
    ax.plot(wavelengths, NInt*spectrum_mid, label='N=%i, Intermediate Mass Stars T=%iK, total'%(NInt,Tmid),c='yellow')#%4500)
    ax.plot(wavelengths, NInt*spectrum_mid_cont, ls=':',c='yellow',lw=.5)
    ax.plot(wavelengths, NSunlike*spectrum_sunlike, label='N=%i, Sun-like Stars T=%iK, total'%(NSunlike,Tsunlike),c='orange')#%5500)
    ax.plot(wavelengths, NSunlike*spectrum_sunlike_cont, ls=':',c='orange',lw=.5)
    ax.plot(wavelengths, NLow*spectrum_low, label='N=%i, Low Mass Stars T=%iK, total'%(NLow,Tlow),c='red')#%7000)
    ax.plot(wavelengths, NLow*spectrum_low_cont, ls=':',c='red',lw=.5)
elif cumul:
    redline = NLow*spectrum_low
    orangeline = NSunlike*spectrum_sunlike+redline#NLow*spectrum_low_cont
    yellowline = NInt*spectrum_mid+orangeline#NSunlike*spectrum_sunlike_cont+NLow*spectrum_low_cont
    blueline = NBlue*spectrum_blue+yellowline#NInt*spectrum_mid_cont+NSunlike*spectrum_sunlike_cont+NLow*spectrum_low_cont
    

    ax.plot(wavelengths, blueline, label='N=%i, Blue, High Mass Stars T=%iK, total'%(NBlue,Tblue),c='cyan')
    ax.fill_between(wavelengths, yellowline, blueline, alpha=.5,color='cyan')

    ax.plot(wavelengths, yellowline, label='N=%i, Intermediate Mass Stars T=%iK, total'%(NInt,Tmid),c='yellow')
    ax.fill_between(wavelengths,orangeline,yellowline,alpha=.5,color='yellow')

    ax.plot(wavelengths, orangeline, label='N=%i, Sun-like Stars T=%iK, total'%(NSunlike,Tsunlike),c='orange')
    ax.fill_between(wavelengths,redline,orangeline,alpha=.5,color='orange')

    ax.plot(wavelengths, redline, label='N=%i, Low Mass Stars T=%iK, total'%(NLow,Tlow),c='red')
    ax.fill_between(wavelengths, 0, redline, alpha=.5,color='red')

else:
    ax.plot(wavelengths, spectrum_blue, label='N=%i, Blue, High Mass Stars T=%iK'%(NBlue,Tblue),c='cyan')#%3000)
    ax.plot(wavelengths, spectrum_mid, label='N=%i, Intermediate Mass Stars T=%iK'%(NInt,Tmid),c='yellow')#%4500)
    ax.plot(wavelengths, spectrum_sunlike, label='N=%i, Sun-like Stars T=%iK'%(NSunlike,Tsunlike),c='orange')#%5500)
    ax.plot(wavelengths, spectrum_low, label='N=%i, Low Mass Stars T=%iK'%(NLow,Tlow),c='red')#%7000)
ax.legend(loc=1,fontsize=12,bbox_to_anchor=(1, -0.1))
ax.plot(wavelengths, galsepect*1.01, label='Simulated Galaxy Spectrum', color='grey',lw=2,zorder=1)
if wemis:
    ax.text(6563+15,np.amax(galsepect)*.99,'H-alpha Emission',color='grey')
    ax.text(5007+15,np.amax(galsepect)*.99,'[OIII] Emission',color='grey')
ax.grid()
c1='k';c2='grey'
ax.yaxis.label.set_color(c1)
ax.tick_params(axis='y', colors=c1)
ax2.yaxis.label.set_color(c2)
ax2.tick_params(axis='y', colors=c2)
ax.grid(c=c1,alpha=.4)
ax2.grid(ls=':',c=c2,alpha=.3)
if summed or cumul:
    ax2.set_ylim(ax.get_ylim())
plt.title('Integrated Galaxy Spectrum from Summed Populations',c=c2)
ax.set_xlabel("    Wavelength (Angstroms)",loc='left')

# plt.yscale('log')

# %% 
#TODO make interactive, 
# probably make the cont/wabs so that you can switch between by making the plotted line the parameter and either continuum line or the spectrum the selection that is toggled
# make the cumulative, summed, separate options toggles
# make emission and absorption lines check boxes to add, on by default 
# sliders for the stellar populations
# also need to make the distrbution width physically motivated instead of mock right now
