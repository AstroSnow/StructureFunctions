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
from scipy.fftpack import fftfreq, ifft2, fft2

#Define a filter here
def box_filter(U,fmin,fmax):
    
    return#(Fro)

#file name
fname='/home/bs428/Documents/PIPrl/KHI2Dr1e2_hr/'
#fname='/media/ben/datadrive1/KHIrl/KHI2Dr1e2_2k/'

#Load in some data
ds=PIPpy.pipread(fname,12,vararrin=['ro_p'])

#Something to emphasise the structures since there is a huge density jump
ro=np.log10(ds['ro_p'][6:-7,6:-7])

#It probably makes sense to do the filtering in the frequency domain so Fourier transform the data
Fro=fft2(ro)

FreqCompRows = fftfreq(Fro.shape[0],d=ds['xgrid'][1]-ds['xgrid'][0])
FreqCompCols = fftfreq(Fro.shape[1],d=ds['ygrid'][1]-ds['ygrid'][0])
#rohat=box_filter(ro,10,100)

#Filter the data based on wavelength
FroFilter=Fro[np.where(FreqCompRows > 10)[0],np.where(FreqCompCols > 10)[0]]

plt.pcolormesh(ro)
plt.show()
