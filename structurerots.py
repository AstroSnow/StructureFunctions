#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:00:09 2022

@author: ben
"""
def structure_function(rvx,rvy,x,y,nsamps,ns):
    import numpy as np
    #ns=100

    sfx=np.zeros(ns)
    sfy=np.zeros(ns)

    for i in range(1,int(nsamps)):
        #i2=np.int64(np.rint((np.size(x)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        #j2=np.int64(np.rint((np.size(y)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        i2=np.int64(np.random.uniform(low=1.0,high=np.size(x)-ns,size=1)[0])
        j2=np.int64(np.random.uniform(low=1.0,high=np.size(y)-ns,size=1)[0])
        #print(i2,j2)
        ta=(rvx[j2,i2:i2+ns+1]-rvx[j2,i2])**3.0
        sfx=sfx+ta[1:ns+1]
        ta=(rvy[j2:j2+ns+1,i2]-rvy[j2,i2])**3.0
        sfy=sfy+ta[1:ns+1]
        
    sfx=sfx/nsamps
    sfy=sfy/nsamps
    
    yr=np.linspace(1, ns,num=ns)*(y[1]-y[0])
    xr=np.linspace(1, ns,num=ns)*(x[1]-x[0])
    
    data={}
    data['sfx']=sfx
    data['sfy']=sfy
    data['xr']=xr
    data['yr']=yr
    
    return(data)

def structure_function_full(rvx,rvy,x,y,nsamps):
    import numpy as np
    #ns=100

    nsx=np.size(rvx[0,:])
    nsy=np.size(rvx[:,0])

    sfx=np.zeros(np.int(np.round((nsx-1)/2)))
    sfy=np.zeros(np.int(np.round((nsy-1)/2)))

    for i in range(1,int(nsamps)):
        #i2=np.int64(np.rint((np.size(x)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        #j2=np.int64(np.rint((np.size(y)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        i2=np.int64(np.random.uniform(low=1.0,high=np.size(x),size=1)[0])
        j2=np.int64(np.random.uniform(low=1.0,high=np.size(y),size=1)[0])
        iarr=[i for i in range(1,nsx)]
        iarr=np.array(iarr)+i2
        iarrs=np.where(iarr<nsx,iarr,iarr-nsx)
        jarr=[i for i in range(1,nsy)]
        jarr=np.array(jarr)+j2
        jarrs=np.where(jarr<nsy,jarr,jarr-nsy)
        #print(iarrs)
        #print(i2,j2)
        ta=(rvx[j2,iarrs]-rvx[j2,i2])**3.0
        sfx=sfx+ta[0:np.int(np.round((nsx-1)/2))]
        ta=(rvy[jarrs,i2]-rvy[j2,i2])**3.0
        sfy=sfy+ta[0:np.int(np.round((nsy-1)/2))]
        
    sfx=sfx/nsamps
    sfy=sfy/nsamps
    
    yr=np.linspace(1, nsy,num=np.int(np.round((nsy-1)/2)))*(y[1]-y[0])
    xr=np.linspace(1, nsx,num=np.int(np.round((nsx-1)/2)))*(x[1]-x[0])
    
    data={}
    data['sfx']=sfx
    data['sfy']=sfy
    data['xr']=xr
    data['yr']=yr
    #return(iarrs)
    return(data)

def structure_function_o3(rv,rvx,rvy,x,y,nsamps):
    import numpy as np
    #ns=100

    nsx=np.size(rvx[0,:])
    nsy=np.size(rvx[:,0])

    sfx=np.zeros(np.int(np.round((nsx-1)/2)))
    sfy=np.zeros(np.int(np.round((nsy-1)/2)))

    for i in range(1,int(nsamps)):
        #i2=np.int64(np.rint((np.size(x)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        #j2=np.int64(np.rint((np.size(y)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        i2=np.int64(np.random.uniform(low=1.0,high=np.size(x),size=1)[0])
        j2=np.int64(np.random.uniform(low=1.0,high=np.size(y),size=1)[0])
        iarr=[i for i in range(1,nsx)]
        iarr=np.array(iarr)+i2
        iarrs=np.where(iarr<nsx,iarr,iarr-nsx)
        jarr=[i for i in range(1,nsy)]
        jarr=np.array(jarr)+j2
        jarrs=np.where(jarr<nsy,jarr,jarr-nsy)
        #print(iarrs)
        #print(i2,j2)
#        ta=(rvx[j2,iarrs]-rvx[j2,i2])**3.0
        ta=(rv[j2,iarrs]-rv[j2,i2])**2.0*(rvx[j2,iarrs]-rvx[j2,i2])
        sfx=sfx+ta[0:np.int(np.round((nsx-1)/2))]
        #ta=(rvy[jarrs,i2]-rvy[j2,i2])**3.0
        ta=(rv[jarrs,i2]-rv[j2,i2])**2.0*(rvy[jarrs,i2]-rvy[j2,i2])
        sfy=sfy+ta[0:np.int(np.round((nsy-1)/2))]
        
    sfx=sfx/nsamps
    sfy=sfy/nsamps
    
    yr=np.linspace(1, nsy,num=np.int(np.round((nsy-1)/2)))*(y[1]-y[0])
    xr=np.linspace(1, nsx,num=np.int(np.round((nsx-1)/2)))*(x[1]-x[0])
    
    data={}
    data['sfx']=sfx
    data['sfy']=sfy
    data['xr']=xr
    data['yr']=yr
    #return(iarrs)
    return(data)
    
def structure_function_scalar(rv,x,y,nsamps):
    import numpy as np
    #ns=100

    nsx=np.size(rv[0,:])
    nsy=np.size(rv[:,0])

    sfx=np.zeros(int(np.round((nsx-1)/2)))
    sfy=np.zeros(int(np.round((nsy-1)/2)))

    for i in range(1,int(nsamps)):
        #i2=np.int64(np.rint((np.size(x)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        #j2=np.int64(np.rint((np.size(y)-ns-1)*np.random.uniform(low=0.0,high=1.0,size=1))[0])
        i2=np.int64(np.random.uniform(low=1.0,high=np.size(x),size=1)[0])
        j2=np.int64(np.random.uniform(low=1.0,high=np.size(y),size=1)[0])
        iarr=[i for i in range(1,nsx)]
        iarr=np.array(iarr)+i2
        iarrs=np.where(iarr<nsx,iarr,iarr-nsx)
        #iarrs=iarr[np.where(iarr<nsx-1)]
        jarr=[i for i in range(1,nsy)]
        jarr=np.array(jarr)+j2
        jarrs=np.where(jarr<nsy,jarr,jarr-nsy)
        #jarrs=jarr[np.where(jarr<nsy-1)]
        #print(iarrs)
        #print(i2,j2)
#        ta=(rvx[j2,iarrs]-rvx[j2,i2])**3.0
        ta=(rv[j2,i2]-rv[j2,iarrs])**3.0
        sfx=sfx+ta[0:int(np.round((nsx-1)/2))]
        
    sfx=sfx/nsamps
    #sfy=sfy/nsamps
    
    yr=np.linspace(1, nsy,num=int(np.round((nsy-1)/2)))*(y[1]-y[0])
    xr=np.linspace(1, nsx,num=int(np.round((nsx-1)/2)))*(x[1]-x[0])
    
    data={}
    data['sfx']=sfx
    #data['sfy']=sfy
    data['xr']=xr
    data['yr']=yr
    #return(iarrs)
    return(data)
    
def structure_function_arb(vx,vy,x,y,nsamps):
    #Compute the structure function for an arbitrary length r using the v_parr v^2 definition
    import numpy as np
    #ns=100

    nsx=np.size(rv[0,:])
    nsy=np.size(rv[:,0])

    sf=np.zeros(int(np.round((nsx-1)/2))) 

    for i in range(1,int(nsamps)):
        i2=np.int64(np.random.uniform(low=1.0,high=np.size(x),size=1)[0])
        j2=np.int64(np.random.uniform(low=1.0,high=np.size(y),size=1)[0])
        iarr=[i for i in range(1,nsx)]
        iarr=np.array(iarr)+i2
        iarrs=np.where(iarr<nsx,iarr,iarr-nsx)
        #iarrs=iarr[np.where(iarr<nsx-1)]
        jarr=[i for i in range(1,nsy)]
        jarr=np.array(jarr)+j2
        jarrs=np.where(jarr<nsy,jarr,jarr-nsy)
        print('NEED TO CALCULATE THE r')
        #r=
        #jarrs=jarr[np.where(jarr<nsy-1)]
        #print(iarrs)
        #print(i2,j2)
#        ta=(rvx[j2,iarrs]-rvx[j2,i2])**3.0
        ta=(v2[j2,i2]-v2[jarrs,iarrs])**2.0
        sfx=sfx+ta[0:int(np.round((nsx-1)/2))]
        
    sfx=sfx/nsamps
    #sfy=sfy/nsamps
    
    yr=np.linspace(1, nsy,num=int(np.round((nsy-1)/2)))*(y[1]-y[0])
    xr=np.linspace(1, nsx,num=int(np.round((nsx-1)/2)))*(x[1]-x[0])
    
    data={}
    data['sfx']=sfx
    #data['sfy']=sfy
    data['xr']=xr
    data['yr']=yr
    #return(iarrs)
    return(data)
    
def spectra_calc(ro):
    import numpy as np
    import scipy.stats as stats
    Fro=np.fft.fft2(ro)
    FroA=np.abs(Fro)**2 #Get magnitude

    npix=ro.shape[0]
    npiy=ro.shape[1]
    FreqCompRows = np.fft.fftfreq(Fro.shape[0])*npix
    FreqCompCols = np.fft.fftfreq(Fro.shape[1])*npiy

    FreqMeshX,FreqMeshy=np.meshgrid(FreqCompCols,FreqCompRows)
    Freq=np.sqrt(FreqMeshX**2+FreqMeshy**2)

    knrm = Freq.flatten()
    fourier_amplitudes = FroA.flatten()

    kbins = np.arange(0.1, npix//2+1, 1.)
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                         statistic = "mean",
                                         bins = kbins)
    Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)
    
    return(Abins,kvals)
