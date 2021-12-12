import matplotlib.pyplot as plt
import numpy as np

absxmax = 100
_dx = 0.05
N = int(2*absxmax/_dx)+1
assert(N%2 == 1) #N//2 is center of space
X = np.linspace(-absxmax,absxmax,N)
assert(len(X) == N)
assert(abs(X[1]-X[0]-_dx) < _dx/10000)
a = 2   
iofxm = N//2+int(a/_dx)+90 # the index, in X, defining X[iofxm]=x_m
mhb2=0.0483                
V0 = 83 # inside the box
E = -50 # E > -V0
ksq0 = -mhb2 * E
tolerance = 0.1
def solve_ode(ksq, forward=True):
    sgn = 1 if forward else -1
    dx = X[1]-X[0]; dx *= sgn
    psiend = np.exp(-np.sqrt(ksq)*absxmax)
    PSI = np.zeros(N)
    dPSI = np.zeros(N)
    PSI[0 if forward else N-1] = psiend
    dPSI[0 if forward else N-1] = np.sqrt(ksq)*psiend * sgn
    it = range(1,N)
    if not forward:
        it = it[::-1]
    for i in it:
        li = i-sgn
        d2psi = ksq*PSI[li] if abs(X[i] > a) else (ksq-mhb2*V0)*PSI[li]
        dPSI[i] = dPSI[li]+d2psi*dx 
        PSI[i] = PSI[li]+dPSI[li]*dx
    return (PSI,dPSI)

def err(pml,pmlp,pmr,pmrp):
    return (pmlp/pml-pmrp/pmr)/(pmlp/pml+pmrp/pmr)

def converge_eps(ksq):
    #step 5.
    pml,pmlp = solve_ode(ksq)
    pmr,pmrp = solve_ode(ksq)
    print(iofxm)
    eps0= err(pml[iofxm],pmlp[iofxm],pmr[iofxm],pmrp[iofxm])
    print(pml,pmr,eps0)
    ksqnew = ksq
    return ksqnew, eps0

kk = ksq0
kk2, eps = converge_eps(ksq0)
while abs(eps) > tolerance:
    break #dbg
    kk = kk2
    kk2, eps = converge_eps(kk)
k2fin = kk2

offset = 2*(iofxm-N//2) 
low = N//2 - offset
high = N//2 + offset
PSI, _  = solve_ode(k2fin)
plt.plot(X[low:high],PSI[low:high])
plt.show()
