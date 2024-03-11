#Spectra test plots
import h5py as h5
import sys
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
sys.path.append('..')
from structurerots import spectra_calc,structure_function_scalar

#Get data
hf=h5.File('DataOZ_magpot/DataOZ/DataOZ_s1.h5','r')
hf=h5.File('DataOZ_512/DataOZ_512_s1.h5','r')
vx=hf['tasks/u'][:,0,:,:,0]
vy=hf['tasks/u'][:,1,:,:,0]
bx=hf['tasks/B'][:,0,:,:,0]
by=hf['tasks/B'][:,1,:,:,0]
b2=bx**2+by**2
v2=vx**2+vy**2
try:
    vz=hf['tasks/u'][:,2,:,:,0]
    bz=hf['tasks/B'][:,2,:,:,0]
    b2+=bz**2
    v2+=vz**2
except:
    pass
#print(hf['tasks/B'].keys())
#x=hf['tasks/B/grid_space']
hf.close

x=np.linspace(0,1,np.shape(vx)[1])
y=np.linspace(0,1,np.shape(vx)[0])
#print(x)

Abins,kvals=spectra_calc(v2[99,:,:])
Abinsb,kvalsb=spectra_calc(b2[99,:,:])


plt.loglog(kvals, Abins,'r')
plt.loglog(kvalsb, Abinsb,'b')
plt.xlabel("$k$")
plt.ylabel("$P(k)$")
plt.tight_layout()
plt.show()


data=structure_function_scalar(b2[90,:,:],x,y,1.0e5)
data2=structure_function_scalar(v2[90,:,:],x,y,1.0e5)

plt.plot(data['xr'],data['sfx'],'b')
plt.plot(data2['xr'],data2['sfx'],'r')
plt.show()
