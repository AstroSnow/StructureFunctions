import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

#hf=h5.File('shearTest/shearTest_s1.h5','r')
#hf=h5.File('DataRecTest/DataRecTest_s1.h5','r')
#hf=h5.File('DataRecNoRhoTest/DataRecNoRhoTest_s1.h5','r')
hf=h5.File('DataOZ/DataOZ_s1.h5','r')
#hf=h5.File('DataOZ1k/DataOZ1k_s1.h5','r')

#rho=hf['tasks/density']
v=hf['tasks/u']
print(np.shape(v))
vx=hf['tasks/u'][:,0,:,:]
vy=hf['tasks/u'][:,1,:,:]
bx=hf['tasks/b'][:,0,:,:]
by=hf['tasks/b'][:,1,:,:]
b2=bx**2+by**2
try:
    vz=hf['tasks/u'][:,2,:,:,0]
    bz=hf['tasks/b'][:,2,:,:,0]
    b2+=bz**2
except:
    pass
hf.close

#print(np.shape(rho))
nx=np.size(vx[0,:,0])
ny=np.size(vx[0,0,:])
X,Y=np.mgrid[0:nx,0:ny]
print(np.shape(X))

#plt.pcolormesh(rho[-1,:,:,0].T)
#plt.pcolormesh(vy[-1,:,:].T)
#plt.pcolormesh(b2[0,:,:].T)
t=30
#plt.pcolormesh(b2[t,:,:])
ss=10
#plt.quiver(Y[0:-1:ss,0:-1:ss],X[0:-1:ss,0:-1:ss],vy[t,0:-1:ss,0:-1:ss],vx[t,0:-1:ss,0:-1:ss])
for t in range(0,100):
    plt.pcolormesh(b2[t,:,:])
    ss=10
    plt.quiver(Y[0:-1:ss,0:-1:ss],X[0:-1:ss,0:-1:ss],vy[t,0:-1:ss,0:-1:ss],vx[t,0:-1:ss,0:-1:ss])
    #plt.plot(rho[i,0,:])
    plt.pause(0.1)
#plt.colorbar()
plt.show()
#for f in hf:
    

