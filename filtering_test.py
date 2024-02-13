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

#Define a filter here?

#Load in some data
#ds=PIPpy.pipread('/media/ben/datadrive1/KHIrl/KHI2Dr1e2_2k/',14,vararrin=['ro_p'])
ds=PIPpy.pipread('/media/snow/datadrive2/KHIrldata/KHI2Dr1e2_hr/',14,vararrin=['ro_p'])

#Something to emphasise the structures since there is a huge density jump
ro=np.log10(ds['ro_p'])
#ro=np.vstack((ro,np.flip(ro))) #mirror the boundaries to make it periodic

#It probably makes sense to do the filtering in the frequency domain so Fourier transform the data
Fro=np.fft.fft2(ro)
FroA=np.abs(Fro) #Get magnitude

FreqCompRows = np.fft.fftfreq(Fro.shape[0])
FreqCompCols = np.fft.fftfreq(Fro.shape[1])

FreqMeshX,FreqMeshy=np.meshgrid(FreqCompRows,FreqCompRows)
Freq=np.sqrt(FreqMeshX**2+FreqMeshy**2)

FroFilter=np.copy(Fro)
#FroFilter[np.where(Freq>0.05)]=0
FroFilter=Fro*np.exp(-(Freq-0.05)**2/0.000001)

roFilter=np.fft.ifft2(FroFilter)

plt.pcolormesh(np.abs(roFilter))