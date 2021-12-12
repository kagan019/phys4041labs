import numpy as np
dgesv = np.linalg.solve
import matplotlib.pyplot as plt

def tridiag(N, a):
    A = np.zeros((N+1,N+1))
    for i in range(N+1):
        A[i][i] = 1+2*a
        if i < N:
            A[i+1][i] = -a
            A[i][i+1] = -a
    return A

def ucol_init(N):
    u = np.zeros(N+1)
    u[int(N/2)] = 1
    return u

def main():
    #a = D(dt/(dx)^2)
    a = 0.01
    N = 25 # number of steps in x
    Nt = 25 # number of t steps to show
    A = tridiag(N, a)
    print(A)
    b = ucol_init(N)
    timesteps = [b]
    for n in range(N+1):
        b = dgesv(A,b)
        b[0], b[N] = 0, 0 #boundary conditions
        timesteps += [b]

    for i,xx in enumerate(timesteps[:Nt]):
        plt.plot(range(len(xx)),xx,label="timestep={}".format(i))
    plt.xlabel("discretized x units")
    plt.ylabel("u(t)")
    plt.legend()
    plt.show()
if __name__ == "__main__":
    main()