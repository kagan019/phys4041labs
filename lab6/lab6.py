# dx/dt = p            dy/dt = q
# dp/dt = -x/r**3 dq/dt = -y/r**3
# x(0) = 0.5 y(0) = 0 p(0) = 0 q(0) = 1.63

import numpy as np
import matplotlib.pyplot as plt

lowt, hight = 0,10
N = 1000000 # number of points
tspace = np.linspace(lowt,hight,N)
sys0 = [[0.5,0],[0,1.63]]
dt = tspace[1]-tspace[0]

def nextsys(xsys):
    xs, ys = xsys
    x, p = xs
    y, q = ys
    r3 = (x**2+y**2)**(3/2)
    dpdt = -x/r3
    dqdt = -y/r3
    return [[x+p*dt, p+dpdt*dt], [y+q*dt, q+dqdt*dt]]

def euler():
    Xs = []
    for i,xi in enumerate(tspace):
        Xs += [sys0] if i == 0 else [nextsys(Xs[i-1])]
    return Xs

def main():
    fpoints = euler()
    xx = [[point[0][0]] for point in fpoints]
    yy = [[point[1][0]] for point in fpoints]
    plt.plot(xx,yy, label="forwards Euler")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
