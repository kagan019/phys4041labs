

#pT/pt = -pT/px + (p^2T/px^2+p^2T/py^2)
#use ADI
#simulate until pT/pt=0
#forward time, centered space
#obtain (x,y,T) for all t:
# (0,y,65), (30,y,25), (x,10,25)
# pT/py = 0 at y = 0

# discretize T. i selects column along x and j selects row along y
# x = 0 --> i = 0 --> T[0][j] = 65
# x = 30 --> i = 30/h --> T[30/h][j] = 25
# y = 10 --> j = 10/h --> T[i][10/h] = 25

EPS=0.08 #lower this for less generous tolerance. But, it might take awhile to compute.

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

N = 700 # i iterated
M = 900 # j iterated
dt = 0.05
h = 0.05
width = N*h
height = M*h
alpha = dt/(2*h**2)
beta = dt/(4*h)
a1 = alpha+beta
b1 = 1+2*alpha
c1 = -(beta+alpha)
a2 = -alpha
b2 = 1+2*alpha
c2 = -alpha
L = N*M

def wrapcol(T): # same-x values are adjacent to eachother in the resulting serialization
    v = np.zeros(N*M)
    for i in range(1,N+1):
        for j in range(1,M+1):
            l = (j-1)*N+i
            v[l-1] = T[i-1][j-1]
    return v
def wraprow(T): # same-y values are adjacent to eachother in the resulting serialization
    N = len(T)
    M = len(T[0])
    w = np.zeros(N*M)
    for i in range(1,N+1):
        for j in range(1,M+1):
            m = (i-1)*M+j
            w[m-1] = T[i-1][j-1]
    return w
def unwrapcol(v):
    T = np.zeros((N,M))
    for l in range(1,L+1):
        #I think the notes are incorrect
        #If they were right, there would be N values of j here
        i = 1 + (l-1)%N
        j = 1 + int((l-1)/N)
        T[i-1][j-1] = v[l-1]
    return T
def unwraprow(w):
    T = np.zeros((N,M))
    for m in range(1,L+1):
        #also adjusted here, so there are N, not M, values of i
        i = 1 + int((m-1)/M)
        j = 1 + (m-1)%M
        T[i-1][j-1] = w[m-1]
    return T

def dcol(T):
    v = wrapcol(T)
    diags = [alpha*np.ones(L-1),(1-2*alpha)*np.ones(L),alpha*np.ones(L-1)]
    t = sp.sparse.diags(diags,(-1,0,1))
    d = t.dot(v)
    #impose boundary conditions
    for j in range(M):
        # for T(x=0) = 65
        l = (j)*N+0
        d[l] = 65
        # for T(x=30) = 25
        l2 = (j)*N+int(30/h)
        d[l2] = 25
    for i in range(N):
        # for T(y=10) = 25 
        l = int(10/h)*N+i
        d[l] = 25 
        # for pT/py|y=0 = 0 
        l2 = 1*N+i
        d[l2] = 0
    return d
def drow(T):
    v = wraprow(T)
    diags = [(alpha+beta)*np.ones(L-1),(1-2*alpha)*np.ones(L),(alpha-beta)*np.ones(L-1)]
    t = sp.sparse.diags(diags,(-1,0,1))
    d = t.dot(v)
    #impose boundary conditions
    for j in range(M):
        # for T(x=0) = 65
        l = 0*M+j
        d[l] = 65
        # for T(x=30) = 25
        l2 = int(30/h)*M+j
        d[l2] = 25
    for i in range(N):
        # for T(y=10) = 25 
        l = i*M+int(10/h)
        d[l] = 25 
        # for pT/py|y=0 = 0 
        l2 = i*M+0
        d[l2] = 0
    return d
def firsthalfu(T):
    d = dcol(T)
    diags = [a1*np.ones(L-1),b1*np.ones(L),c1*np.ones(L-1),np.zeros(L-N)]
    #impose boundary conditions
    def seteqncoeff(l,dg,lt,rt):
        l = N*j+0
        diags[1][l] = dg
        if l != 0:
            diags[0][l-1] = lt
        if l != L-1:
            diags[2][l] = rt
    for j in range(M):
        # for T(x=0) = 65
        l = N*j+0
        seteqncoeff(l,1,0,0)
        # for T(x=30) = 25
        l2 = N*j+int(30/h)
        seteqncoeff(l,1,0,0)
    for i in range(N):
        # for T(y=10) = 25
        l = N*int(10/h)+i
        seteqncoeff(l,1,0,0)
        # for pT/py|y=0 = 0
        l2 = N*1+i
        seteqncoeff(l2,1,0,0)
        diags[3][l2-N] = -1
    tridiag = sp.sparse.diags(diags, (-1,0,1,-N),format="csc")
    return sp.sparse.linalg.spsolve(tridiag, d)
def secondhalfu(d):
    d = drow(T)
    diags = [a2*np.ones(L-1),b2*np.ones(L),c2*np.ones(L-1)]
    #impose boundary conditions
    def seteqncoeff(l,dg,lt,rt):
        l = 0*M+j
        diags[1][l] = dg
        if l != 0:
            diags[0][l-1] = lt
        if l != L-1:
            diags[2][l] = rt
    for j in range(M):
        # for T(x=0) = 65
        l = 0*M+j
        seteqncoeff(l,1,0,0)
        # for T(x=30) = 25
        l2 = int(30/h)*M+j
        seteqncoeff(l,1,0,0)
    for i in range(N):
        # for T(y=10) = 25
        l = i*M+int(10/h)
        seteqncoeff(l,1,0,0)
        # for pT/py|y=0 = 0
        l2 = i*M+0
        seteqncoeff(l2,1,-1,0)
    tridiag = sp.sparse.diags(diags, (-1,0,1),format="csc")
    return sp.sparse.linalg.spsolve(tridiag, d)

Ts = []
T = np.random.rand(N,M)*100
Tfull = np.array(T)
three =0
while three<3 or abs(sum(wrapcol(Tfull) - wrapcol(T)) / (N*M)) > EPS:
    #ensure the animations evolves past initial noise
    three += 1
    three = 3 if three > 3 else three
    # this term is observed to strictly decrease, after the first 
    # couple iterations, for ~10000 iterations. However, marginal
    # decreases gradually slow as the simulation goes on.
    print(abs(sum(wrapcol(Tfull) - wrapcol(T)) / (N*M))) 
    # confirm boundary conditions.
    # we actually implement them elsewhere though. so except for the first
    # iteration this part is redundant
    for m in range(M):
        T[0][m] = 65
        T[int(30/h)][m] = 25
        T[0][m] = T[1][m]
    for n in range(N):
        T[n][int(10/h)] = 25
        T[n][1] = T[n][0]
    Ts += [T]
    T = Tfull

    uhalf = firsthalfu(T)
    Thalf = unwrapcol(uhalf)
    ufull = secondhalfu(drow(Thalf))
    Tfull = unwraprow(ufull)

#plot
fig,ax = plt.subplots()
def animate(i):
    ax.pcolormesh(Ts[i],cmap=plt.cm.get_cmap('jet',100), rasterized=True)

anim = FuncAnimation(fig, animate, interval=1000 , repeat=True,frames=len(Ts))

plt.xlabel("xpos")
plt.ylabel("ypos")
plt.draw()
plt.show()
