import numpy as np
import sys

file = sys.argv[1]
out = file.replace(".dat",".out")

nbin = len([1 for line in open(file)])

N_meas = 4
data = np.loadtxt(file).reshape((nbin,-1,N_meas+1))

mean = np.nanmean(data,axis=0)
error = np.nanstd(data,axis=0)/np.sqrt(nbin)

S = mean[:,0]
mean = mean[:,1:]
error = error[:,1:]

N_row = mean.shape[0]

proc_data = np.zeros((N_row,2*N_meas+1),dtype=np.float64)

proc_data[:,0] = S
proc_data[:,1::2] = mean
proc_data[:,2::2] = error


np.savetxt(out,proc_data,fmt="%20.10e")

