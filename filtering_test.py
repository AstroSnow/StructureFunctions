#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for applying filtering to some data to extract scales
Created on Thu Feb  8 09:41:28 2024
@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods as PIPpy
import scipy.stats as stats

#Define a filter here?

#Load in some data
#ds=PIPpy.pipread('/media/ben/datadrive1/KHIrl/KHI2Dr1e2_2k/',14,vararrin=['ro_p']); ro=np.log10(ds['ro_p'])
ds=PIPpy.pipread('/media/snow/datadrive2/KHIrldata/KHI2Dr1e2_hr/',6,vararrin=['ro_p']); ro=np.log10(ds['ro_p'][6:-7,6:-7]);ro=np.vstack((ro,np.flip(ro)))
#ds=PIPpy.pipread('/media/snow/store1/reconnection_data/mhd_rad_8k/',40,vararrin=['ro_p']); ro=ds['ro_p'][1:-1,1:-1]

#Something to emphasise the structures since there is a huge density jump
#ro=np.vstack((ro,np.flip(ro))) #mirror the boundaries to make it periodic

#It probably makes sense to do the filtering in the frequency domain so Fourier transform the data
Fro=np.fft.fft2(ro)
FroA=np.abs(Fro)**2 #Get magnitude

npix=ro.shape[0]
npiy=ro.shape[0]
FreqCompRows = np.fft.fftfreq(Fro.shape[0])*npix
FreqCompCols = np.fft.fftfreq(Fro.shape[1])*npiy

FreqMeshX,FreqMeshy=np.meshgrid(FreqCompCols,FreqCompRows)
Freq=np.sqrt(FreqMeshX**2+FreqMeshy**2)

FroFilter=np.copy(Fro)
#FroFilter[np.where(Freq>10.0)]=0
#FroFilter=Fro*np.exp(-(Freq-10.0)**2/100.0)
FroFilter=Fro*(0.5*(np.tanh((Freq-50)/1.0)+1)-0.5*(1+np.tanh((Freq-200)/10.0)))


knrm = Freq.flatten()
fourier_amplitudes = FroA.flatten()
fourier_amplitudes_filter = np.abs(FroFilter**2).flatten()

kbins = np.arange(0.1, npix//2+1, 1.)
kvals = 0.5 * (kbins[1:] + kbins[:-1])
Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                     statistic = "mean",
                                     bins = kbins)
Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

AbinsFilter, _, _ = stats.binned_statistic(knrm, fourier_amplitudes_filter,
                                     statistic = "mean",
                                     bins = kbins)
AbinsFilter *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

fig,ax=plt.subplots()
ax.loglog(kvals,Abins,'k')
ax.loglog(kvals,AbinsFilter,'r')
ax.set_ylim(np.min(Abins),np.max(Abins))

roFilter=np.fft.ifft2(FroFilter)

fig,ax=plt.subplots()
#plt.pcolormesh(np.abs(ro))
ax.pcolormesh(np.abs(roFilter))