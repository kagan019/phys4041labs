import numpy as np
from timeit import default_timer
import matplotlib.pyplot as plt


def rand_matrix(n):
    return np.array([
        [np.random.rand() for _ in range(n)] 
        for _ in range(n)]
    ) + np.identity(n)

def rand_column_vec(n):
    return np.array(
        [np.random.rand() for _ in range(n)]
    ).T

def benchmark(x):
    size= int(x)
    print(size)
    timsum = 0
    iters = 50
    for _ in range(iters):
        start = default_timer()
        np.linalg.solve(rand_matrix(size), rand_column_vec(size))
        end = default_timer() 
        timsum += end-start
    return timsum/iters

def main():
    seed = 82378423
    np.random.seed(seed)
    sizes = np.array(range(1,300,10))
    plt.plot(sizes, np.array([benchmark(z) for z in sizes]))
    plt.xlabel("matrix size, logarithmic")
    plt.ylabel("time")
    plt.ylabel("avg time for 50 solutions")
    plt.xscale('log')
    plt.show()

if __name__ == "__main__":
    main()