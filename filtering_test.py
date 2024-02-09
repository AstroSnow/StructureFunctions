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
ds=PIPpy.pipread('/media/ben/datadrive1/KHIrl/KHI2Dr1e2_2k/',14,vararrin=['ro_p'])

#Something to emphasise the structures since there is a huge density jump
ro=np.log10(ds['ro_p'])

#It probably makes sense to do the filtering in the frequency domain so Fourier transform the data
Fro=np.fft.fft2(ro)
