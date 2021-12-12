import numpy as np
import numpy.linalg as la

A = np.array([[3,2],[2,6]])
b = np.array([[2],[-8]])

def main():
    print(la.solve(A,b))

if __name__ == "__main__":
    main()
