# I'll start with 2. 
# psi(x, t=0) is an eigenstate of the hamiltonian (inside well) iff
# H(psi) = E*psi
# iff (-(hbar^2/2m)(d^2/dx^2)+V(x,t=0))psi = k0^2/(2m)*psi
#   hbar^2/2m = 1
#   V(x) = 0
#   d^2/dx^2 psi = (-1/2*(2x-10)+ik0)^2 * psi 
#       = (x^2 + 10x + 2ixk0 + 25 + (10i-1)k0) * psi
# iff d^2/dx^2 psi = k0^2/(2m) * psi iff false
# so, it is not an eigenstate

#1. 
import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt

k0 = 17*np.pi
sigma0 = 0.5
dx = 0.02; dt = dx**2/2
h = dt/(2*dx**2)
Nx = int(15/dx)+1; Nt = int(1/dt)+1
X = [np.linspace(0,15,Nx) for _ in range(Nt)]
def Psi0(x):
    return np.exp(-1/2*((x-5)**2+sigma0**2)+k0*x*1j)
def RePsi0(x):
    return np.real(Psi0(x))
def ImPsi0(x):
    return np.imag(Psi0(x))
PSI = [[Psi0(x) for x in X[0]]]
PSI[0][0] = 0; PSI[0][Nx-1] = 0

# know: PSI[t][0] = PSI[t][Nx] = 0 for all t

# du/dt = Dd^2u/dx^2 for D = -1/i = i, u = psi
# dpsi/dt = id^2psi/dx^2
#   dpsi/dt = Re(dpsi/dt) + Im(dpsi/dt)i
#   d^2psi/dx^2 = Re(d^2psi/dx^2) + Im(d^2psi/dx^2)i
# Re(dpsi/dt) + Im(dpsi/dt)i = i(Re(d^2psi/dx^2) + Im(d^2psi/dx^2)i)
# we have a system of 2 eqs to solve
# Re(dpsi/dt) = -Im(d^2psi/dx^2) 
# Im(dpsi/dt) = Re(d^2psi/dx^2) 

# Crank-Nicholson:
# u[n+1][j] - u[n][j]
#   = Ddt/(2dx^2)[
#     (u[n+1][j+1]-2u[n+1][j]+u[n+1][j-1])
#     +(u[n][j+1]-2u[n][j]+u[n][j-1])
#   ]
# find u[n+1]
#   h = dt/(2dx^2)
#   U[n] = the column vector { u[n][j] } for all points in space j
#   D = i
# -ihu[n+1][j+1] +(1+2ih)u[n+1][j] -ihu[n+1][j-1] 
#   = ihu[n][j+1]+(1-2ih)u[n][j]+ihu[n][j-1]
# tridiag(-ih,1+2ih,-ih)U[n+1] = tridiag(ih,1-2ih,ih)U[n]
#    Q                                P

def init_tridiagQ():
    diags = [-1j*h*np.ones(Nx-3),(1+2j*h)*np.ones(Nx-2),-1j*h * np.ones(Nx-3)]
    return sp.sparse.diags(diags,(-1,0,1), format="csr")

def init_tridiagP():
    diags = [1j*h*np.ones(Nx-3),(1-2j*h)*np.ones(Nx-2),1j*h * np.ones(Nx-3)]
    return sp.sparse.diags(diags,(-1,0,1), format="csr")
tQ = init_tridiagQ()
tP = init_tridiagP()
def cn(Un):
    # returns U[n+1], excluding the endpoints which are 0
    b = np.matmul(tP.toarray(),np.array(Un))
    return sp.sparse.linalg.spsolve(tQ,b)

print("to {}".format(len(X)))
for i,x in enumerate(X):
    if i % 1000 == 0:
        print(i)
    nextU = cn(PSI[i][1:-1])
    nextU =  np.concatenate(([0],nextU,[0]))
    PSI += [nextU]



# 3. plot ||psi||^2 vs x
numsecs = int(dt*(Nt-1)+1)
#timestep_indexes = np.array([int(a) for a in np.array(range(numsecs))/dt])
timestep_indexes = np.array([0, int(0.1/dt),int(0.2/dt),int(0.3/dt)])
for t in timestep_indexes:
    plt.plot(X[int(t)],[psi*np.conj(psi) for psi in PSI[int(t)]],label="seconds = {}".format(t*dt))
plt.legend()
plt.show()