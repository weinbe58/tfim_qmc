import numpy as np
import sys

file = sys.argv[1]
out = file.replace(".dat",".out")

nbin = len([1 for line in open(file)])


data = np.loadtxt(file).reshape((nbin,-1))
S_list = np.unique(data[:,0])
Ss = data[:,0]

mean = np.vstack((data[Ss==S,1:].mean(axis=0)) for S in S_list)
error  = np.vstack((data[Ss==S,1:].std(axis=0) )/np.sqrt((Ss==S).sum()) for S in S_list)

proc_shape = S_list.shape
proc_shape += (2*mean.shape[1]+1,)
proc_data = np.zeros(proc_shape,dtype=mean.dtype)


proc_data[:,0] = S_list
proc_data[:,1::2] = mean
proc_data[:,2::2] = error

np.savetxt(out,proc_data,fmt="%30.15f")



