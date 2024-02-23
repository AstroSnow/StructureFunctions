import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

hf=h5.File('test/test_s1.h5','r')

rho=hf['tasks/u'][90,0,:,:]
hf.close

plt.pcolormesh(rho.T)
plt.show()
#for f in hf:
    

