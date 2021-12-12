# dx/dt = x - xy
# dy/dt = -y + xy
# x(0) = 0.1 y(0) = 0.1
# I choose to rewrite this in a more concise way
#   deltaX = X[i+1] - X[i] = [[1, -x`], [y`, -1]]X` * deltat 
#       for x`,y` (and thus X`) at some point in [X[i],X[i+1]]
# So we can choose to approximate deltaX with X` = X[i] 
# (forwards Euler) or X` = X[i+1] (backwards Euler)
# Call the matrix in the rhs A(x`,y`) = A`.

# I choose the forwards Euler for simplicity:
#    X[i+1] = (I + A[i])X[i] * deltat
# let V = (I + A`) = [[2,-x`],[y`,0]]

import numpy as np
import matplotlib.pyplot as plt

lowt, hight = 0,10
N = 100000 # number of points
tspace = np.linspace(lowt,hight,N)
xinit, yinit = 0.1,0.1
dt = tspace[1]-tspace[0]

def V(xt,yt):
    return [[2,-xt], [yt, 0]]

def euler(ts, step_mtx):
    
    def next(xy): 
        mtx = step_mtx(xy[0],xy[1])
        dot = [
            xy[0]*mtx[0][0] + xy[1]*mtx[0][1],
            xy[0]*mtx[1][0] + xy[1]*mtx[1][1]
        ]
        return [dot[0]*dt,dot[1]*dt]
    Xs = []
    for i,xi in enumerate(ts):
        Xs += [[xinit, yinit]] if i == 0 else [next(Xs[i-1])]
    return Xs

def main():
    tspace = np.linspace(lowt,hight,N)
    fpoints = euler(tspace,V)
    xx,yy = [xy[0] for xy in fpoints], [xy[1] for xy in fpoints]
    plt.plot(xx,yy, label="forwards Euler")
    plt.legend()
    plt.show() #the graph is startlingly linear. 
    # This result describes both fish populations starting in equal ratio, and
    # consuming at such pace so that they maintain equal ratio as they steadily
    # decline to mutual extinction.

if __name__ == "__main__":
    main()
