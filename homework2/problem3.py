import numpy as np
from scipy import interpolate


def main():
    query = 3.25
    x = np.array([1.0,1.1,1.2, 1.3,1.4,1.5])
    y = np.array([4.7,3.5,3.0, 2.7,3.0,2.4])-query
    tck = interpolate.splrep(x,y,k=3,s=0)
    solution = interpolate.sproot(tck)[0]
    print("p={} is interpolated to v={}".format(query, solution))

if __name__ == "__main__":
    main()