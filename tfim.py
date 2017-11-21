from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian 
import numpy as np



L = 4
beta=1.0

basis = spin_basis_1d(L)

J = [[-1.0,i,(i+1)%L] for i in range(L)]
h = [[-1.0,i] for i in range(L)]
m = [[1.0/L,i] for i in range(L)]

dynamics = [["zz",J,lambda x:x,()],["x",h,lambda x:1-x,()]]
H = hamiltonian([],dynamics,basis=basis,dtype=np.float64)
M = hamiltonian([["z",m]],[],basis=basis,dtype=np.float64)
M2 = M**2

for t in np.linspace(0,1,11):
	E,V = H.eigh(time=t)

	p = np.exp(-beta*(E-E[0]))
	p /= p.sum()

	M2_diag = M2.expt_value(V,enforce_pure=True)

	print t,(p*M2_diag).sum()

	

