from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian 
import numpy as np
import scipy.sparse as sp
import sys


def f1(t,t_f,S_f):
	return S_f*(t/t_f)

def f2(t,t_f,S_f):
	return 1-S_f*(t/t_f)

L = int(sys.argv[1])
t_f = L*float(sys.argv[2])
S_f = 1.0
mm = int(t_f*L)


basis = spin_basis_1d(L,pauli=True)

JJ = -1.0
hh = -1.0

J = [[JJ,i,(i+1)%L] for i in range(L)]
JI = [[-np.abs(JJ),i] for i in range(L)]
h = [[hh,i] for i in range(L)]
hI = [[-np.abs(hh),i] for i in range(L)]
m = [[1.0,0]]

dynamics1 = [["zz",J,f1,(t_f,S_f)],["I",JI,f1,(t_f,S_f)]]
dynamics2 = [["x",h,f2,(t_f,S_f)],["I",hI,f1,(t_f,S_f)]]
Hzz = hamiltonian([],dynamics1,basis=basis,dtype=np.float64,check_symm=False,check_herm=False)
Hx = hamiltonian([],dynamics2,basis=basis,dtype=np.float64,check_symm=False,check_herm=False)
M = hamiltonian([["z",m]],[],basis=basis,dtype=np.float64,check_symm=False,check_herm=False)
# M2 = M**2


# psi0 = np.zeros((Hzz.Ns,),dtype=np.float64)
# psi0[0] = 1.0

E0,psi0 = Hx.eigsh(k=1,which="SA")
psi0 = psi0.ravel()

# n=0
# states = [np.random.permutation([0,1]) for i in range(L-n)]
# states.extend((np.array([1,1])/np.sqrt(2) for i in range(n)))
# states = np.array(states)
# print np.random.permutation(states)
# psi0 = reduce(np.kron,states[:])

# psi0 = np.zeros(basis.Ns)
# psi0[0] = 1.0

psileft = M.dot(psi0)
print np.abs(M.matrix_ele(psileft,psi0))
# psi0 = np.identity(basis.Ns)
H=(Hzz+Hx+0.0*M)
for i in range(mm):
	t = t_f*i/mm
	psi0 = H.dot(psi0,time=t)
	norms = np.linalg.norm(psi0,axis=0)
	psi0 = (psi0.T/norms).T

	print np.abs(M.matrix_ele(psileft,psi0))


print L,1.0/t_f,np.abs(M.matrix_ele(psileft,psi0))

