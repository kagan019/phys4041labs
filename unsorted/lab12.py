import numpy as np

N = 100
S = 10000
bound = 5
dx = 0.1

def E(x): #dimensionless energy, without vel
    return (x-1)**2*(x+1)**2
def boltz(energy): #if energy is unitless
    return np.exp(-energy)

xs = []
for i in range(N):
    xi = np.random.choice(range(bound))
    for k in range(S):
        #take data on which state ptcl is in...
        xs += [xi*dx]

        V = E(xi*dx)
        r = np.random.uniform()
        xip = None
        if r < 0.5:
            xip = xi-1
        else:
            xip = xi+1
        Vp = E(xip*dx)

        #metropolis
        r2 = np.random.uniform()
        if Vp <= V:
            xi = xip
        elif r2 < boltz(Vp-V):
            xi = xip

print(np.average(xs))