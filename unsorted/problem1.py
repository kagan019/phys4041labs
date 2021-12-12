import matplotlib.pyplot as plt
import numpy as np

absxmax = 15 #not so far as to run into floating point issues
_dx = 0.05
N = int(2*absxmax/_dx)+1
assert(N%2 == 1) #N//2 is center of space
X = np.linspace(-absxmax,absxmax,N)
assert(len(X) == N)
assert(abs(X[1]-X[0]-_dx) < _dx/10000)
a = 2   
iofxm = N//2+int(a/_dx)+60 # the index, in X, defining X[iofxm]=x_m
mhb2=0.0483                
V0 = 83 # inside the box
E = -20 # E > -V0
ksq0 = -mhb2 * E
tolerance = 0.1
def solve_ode(ksq, forward=True):
    sgn = 1 if forward else -1
    dx = X[1]-X[0]; dx *= sgn
    psiend = np.exp(-np.sqrt(ksq)*absxmax)
    PSI = np.zeros(N)
    dPSI = np.zeros(N)
    PSI[0 if forward else N-1] = psiend
    dPSI[0 if forward else N-1] = -np.sqrt(ksq)*psiend*sgn
    it = list(range(N))
    if not forward:
        it = it[::-1]
    for i in it[1:]:
        li = i-sgn
        d2psi = ksq*PSI[li] if abs(X[i]) > a else (ksq-mhb2*V0)*PSI[li]
        dPSI[i] = dPSI[li]+d2psi*dx 
        PSI[i] = PSI[li]+dPSI[li]*dx
    return (PSI,dPSI)

def err(pml,pmlp,pmr,pmrp):
    return (pmlp/pml-pmrp/pmr)/(pmlp/pml+pmrp/pmr)

def converge_eps(ksq,lstksq):
    #step 5.
    #strategy: generic gradient descent to find the 0 of the error function
    pml,pmlp = solve_ode(ksq)
    pmr,pmrp = solve_ode(ksq,forward=False)
    eps1 = err(pml[iofxm],pmlp[iofxm],pmr[iofxm],pmrp[iofxm])
    pml,pmlp = solve_ode(lstksq)
    pmr,pmrp = solve_ode(lstksq,forward=False)
    eps0 = err(pml[iofxm],pmlp[iofxm],pmr[iofxm],pmrp[iofxm])
    derrdksq = (eps1-eps0)/(ksq-lstksq)
    ksqnew = None
    if np.sign(ksq) == np.sign(derrdksq):
        ksqnew = ksq - 0.0003*np.sign(ksq)
    else:
        ksqnew = ksq + 0.0003 * np.sign(ksq)
    return ksqnew, ksq, eps1

kk = ksq0
kk2 = kk+0.0003
eps = 1000
while abs(eps) > tolerance:
    kk2, kk, eps = converge_eps(kk2,kk)    
k2fin = kk2
print(k2fin) #0.9666
PSIleft, _  = solve_ode(k2fin)
PSIright, _2  = solve_ode(k2fin,forward=False)
plt.plot(X,PSIleft, label="from left")
plt.plot(X, PSIright, label="from right")
plt.legend()
plt.show() #the plot is symmetrical. but boring; I expected to see more oscillations?
# I think part of the issue is the numerical conditioning of the problem; floating point wise,
# psiend is prone to spontaneous changes in sign if it is too small