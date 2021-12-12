#dT/dx dy/dx + Td^2y/dx^2 = pd^2y/dt^2
# -z = dT/dx y + T dy/dx
# s = p(x)dy/dt
#ds/dt = -dz/dx
#dq = dt/(2*dx)
#s[n+1][j] = 
#    1/2*(s[n][j+1]-s[n][j-1])- dq(z[n][j+1]-z[n][j-1])
# dT/dx y + T dy/dx + z = 0
# dr/dt = -d/dx [ F(r,s) ]
# ds/dt = -d/dx [ G(r,s) ]

import numpy as np

dt = 0.01
dx = 0.01
N = int(1 / dx)
p0 = 3
T0 = 3
alpha = 5
def p(x):
    return p0*np.exp(alpha*x)
def T(x):
    return T0* np.exp(alpha*x)
def dTdx(x):
    return alpha*T0*np.exp(alpha*x)
def d2ydt2(x, dx, d2ydx2, dydx):
    return (dTdx(x)*dydx+T(x)*d2ydx2)/p(x)

X = [0]
Y = [0]
TT = [0] 
for i in range(N):
    pass
X += [1]
Y += [0]