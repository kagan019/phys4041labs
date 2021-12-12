#pT/pt = -pT/px + (p^2T/px^2+p^2T/py^2)
#use ADI
#simulate until pT/pt=0
#forward time, centered space
#obtain (x,y,T) for all t:
# (0,y,65), (30,y,25), (x,10,25)
# pT/py = 0 at y = 0

# x = 0 --> i = 0 --> T[0][y] = 65
# x = 30 --> i = 30/h --> T[30/h][y] = 25
# y = 10 --> j = 10/h --> T[x][10/h] = 25


import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg

def wrapcol(T):
    N = len(T)
    M = len(T[0])
    v = np.zeros(N*M)
    for i in range(1,N+1):
        for j in range(1,M+1):
            l = (j-1)*N+i
            v[l-1] = T[i-1][j-1]
    return v
def wraprow(T):
    N = len(T)
    M = len(T[0])
    w = np.zeros(N*M)
    for i in range(1,N+1):
        for j in range(1,M+1):
            m = (i-1)*M+j
            w[m-1] = T[i-1][j-1]
    return w
def unwrapcol(v, N):
    L = len(v)
    T = np.zeros(N,L/N)
    for l in range(1,L+1):
        i = 1 + int((l-1)/N)
        j = 1 + (l-1)%N
        T[i-1][j-1] = v[l-1]
    return T
def unwraprow(w, M):
    MM = len(w)
    T = np.zeros(MM/M)
    for m in range(1,MM+1):
        i = 1 + (m-1)%MM
        j = 1 + ((m-1)/MM)
        T[i-1][j-1] = w[m-1]
    return T

dt = 0.05
h = 0.05
alpha = dt/(2*h**2)
beta = dt/(4*h)
a = alpha+beta
b = 1+2*alpha
c = beta-alpha
def dcol(T):
    v = wrapcol(T)
    diags = [alpha,1-2*alpha,alpha]
    t = sp.sparse.diags(diags,(-1,0,1))
    return np.matmul(t,v)
def drow(T):
    w = wraprow(T)
    diags = [alpha,1-2*alpha,alpha]
    t = sp.sparse.diags(diags,(-1,0,1))
    return np.matmul(t,w)
def secondhalfu(d):
    d = [a,b,c]
    t = sp.diags(d, (-1,0,1))
    t[0] = [1, 0, 0, 0]
    for all i t[i][10/h] = 1
    return sp.spsolve(t, d)
def firsthalfu(d):
    d = [a,b,c]
    t = sp.diags(d, (-1,0,1))
    return sp.spsolve(t, d)


T= np.array([[65]])
Tfull = np.array(T)
while sum(wrapcol(Tfull-T)) >
    uhalf = firsthalfu(dcol(T))
    Thalf = unwrapcol(uhalf)
    ufull = secondhalfu(drow(Thalf))
    Tfull = unwraprow(ufull)

