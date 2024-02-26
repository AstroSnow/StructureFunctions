import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

#hf=h5.File('shearTest/shearTest_s1.h5','r')
hf=h5.File('DataRecTest/DataRecTest_s1.h5','r')

rho=hf['tasks/density']
vx=hf['tasks/u'][0,:,:,:,0]
vy=hf['tasks/u'][1,:,:,:,0]
bx=hf['tasks/b'][0,:,:,:,0]
by=hf['tasks/b'][1,:,:,:,0]
bz=hf['tasks/b'][2,:,:,:,0]
b2=bx**2+by**2+bz**2
hf.close

print(np.shape(rho))

#plt.pcolormesh(rho[-1,:,:,0].T)
#plt.pcolormesh(vx[-1,:,:].T)
plt.pcolormesh(b2[0,:,:].T)
#for i in range(0,10):
#	plt.plot(rho[i,0,:])
plt.show()
#for f in hf:
    

