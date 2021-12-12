import numpy as np
from numpy import linalg as la

def u(i,x):
    qi = i*np.pi
    return np.sin(qi*x)
def uHu(j,i,x):
    qi = i*np.pi
    return u(j,x) * qi*qi * np.sin(qi*x) + x**2 * u(i,x) / 2
def uu(j,i,x):
    return u(j,x)*u(i,x)


def integ01(fofx):
    #simpson's rule
    samplepoints = 101
    lsp = np.linspace(0,1,samplepoints)
    h = lsp[1]-lsp[0]
    intgl = 0
    alt = 0
    for xi in lsp:
        coef = 2*(alt+1)
        alt += 1; alt %= 2
        intgl += coef * fofx(xi)
    intgl -= fofx(0) + fofx(1)
    intgl *= h/3
    return intgl

def gen_H(n):
    H = np.zeros((n,n))
    S = np.zeros((n,n))
    for j in range(n):
        for i in range(n):
            g = lambda x: uu(j+1,i+1,x)
            f = lambda x: uHu(j+1,i+1,x)
            H[j][i] = integ01(f)
            S[i][j] = integ01(g)
    return H,S

def energies(num_eigenfuncs):
    H,S = gen_H(num_eigenfuncs)
    w,v = la.eig(np.matmul(la.inv(S), H)) 
    return w

def main():
    np.set_printoptions(precision=2, suppress=False)
    levels = energies(32)
    levels = np.array(sorted(levels))
    print(levels)

    #let me show you this
    E1 = levels[0]
    E = []
    for n in range(1, len(levels)):
        E += [E1 * n*n]
    E = np.array(levels)
    
    offset = []
    for i, j in zip(E,levels):
        offset += [i-j]
    print(np.array(offset)) # this is a good sign


if __name__ == "__main__":
    main()

# (H-E)a = 0
# det(H-E) = 0 