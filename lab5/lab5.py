import numpy as np
import scipy.sparse
from scipy.sparse.linalg import eigs

def p(x):
    return (1+x**2)

def discretize(p, x0, xmax, n):
    X = np.linspace(x0,xmax,n-1)
    xi = -1000
    h = X[1] - X[0]
    Y = np.array([])
    M = [[] for _ in n]
    for i in range(1,n):
        y = np.zeros(n)
        Y += y[i-1] + p(X[i]) * 0
    for i in range(n):
        
        
        for j, xipo in enumerate(X):
            if i == 0:
                continue
            yp[]
            dy1 = y[j+1] - y[j]
            dy2 = y[j] - y[j-1]
            M[i][j-1] = (p(X[i] + h/2)*(dy1/h) - p(x-h/2)*(dy2,h))/h
        

    return M
def main():
    d = discretize(p,-1,1,100)
    e = eigs(d)
    print(min(e))

if __name__ == "__main__":
    main()