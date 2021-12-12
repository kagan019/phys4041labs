import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.tan(x)
def fprime(x):
    b = 1/np.cos(x)
    return b*b

def forward_difference_err(xi,xpdx):
    #lowest order approx
    # (f(x+deltax) - f(x))/deltax
    numer_val = (f(xpdx)-f(xi)) / (xpdx-xi)
    analytic_val = fprime(xi)
    return abs((analytic_val-numer_val))/analytic_val

def central_difference_err(xi,xpdx):
    #lowest order approx
    # (f(x+deltax/2) - f(x-deltax/2))/deltax
    dx = (xpdx-xi)/2
    numer_val = (f(xi+dx)-f(xi-dx)) / (2*dx)
    analytic_val = fprime(xi)
    return abs(analytic_val-numer_val)/analytic_val

def benchmark(n):
    low, high = (-2,2)
    ivl = lambda l, h: np.linspace( l, h,  int(n*(h-l)/(high-low)) ) 
    lsps = [ivl(low, -np.pi/2-0.1), ivl(-np.pi/2 + 0.1, np.pi/2 - 0.1), ivl(np.pi/2+0.1, high)]
    maxfderr = -1
    maxcderr = -1
    for lsp in lsps:
        xi = -1000
        for xpdx in lsp:
            if xi < low: #skip first
                xi = xpdx
                continue
            maxfderr = max(forward_difference_err(xi,xpdx),maxfderr)
            maxcderr = max(central_difference_err(xi,xpdx),maxcderr)

            xi = xpdx

    return maxfderr,maxcderr

def main():
    nn = np.array(range(100,1000))
    fderrs = []
    cderrs = [] 
    for n in nn:
        print(n)
        fde,cde = benchmark(n)
        fderrs += [fde]
        cderrs += [cde]
    plt.plot(nn, fderrs, label="forward difference")
    plt.plot(nn,cderrs, label="central difference")
    plt.ylabel("error, as ratio")
    plt.xlabel("# sampled points")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.show() #shows that the fd converges linearly while cd converges quadratically



if __name__ == "__main__":
    main()
