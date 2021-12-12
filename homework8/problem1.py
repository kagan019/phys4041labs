import numpy as np


def estimateI(N):
    M = 1000
    def I():
        #monte carlo method
        #The integrand (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)**2 describes an 11 dimensional
        #volume, with the expression as the 11th. So the chance a random point (11 
        # random uniform dimensions) lies within that volume (ie, the 11th 
        # dimension is less than or equal to the integrand) is the ratio between the 
        # total volume of space and the volume described by the integrand. 
        integrand = lambda x: (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9])**2
        ten_dimensional_pts = np.random.uniform(size=(N,10))
        #integrand can be as large as 9**2 = 81
        #that's the height of the 11 dimensional hyperrectangle
        vol = 81
        random_pts = np.random.uniform(size=N) * vol
        inside = 0
        for i in range(N):
            pt = ten_dimensional_pts[i]
            if integrand(pt) >= random_pts[i]:
                inside += 1
        return inside*vol/N

    Isum = 0
    for _ in range(M):
        Isum += I()
    return Isum/M
    

N = 8
finalN = 32768
iters = int(np.log2(finalN))-int(np.log2(N))+1
Ns = np.zeros(iters)
Is = np.zeros(iters)
i = 0
while N <= 32768:
    Is[i] = estimateI(N)
    Ns[i] = N
    i += 1
    N *= 2

expected_val = 155/6
print(Is)
import matplotlib.pyplot as plt

plt.plot(1/np.sqrt(Ns),[abs(x-expected_val) for x in Is])
plt.xlabel("1/sqrt(N)")
plt.ylabel("I(N) deviation from expected")
plt.show() # this looks good; larger N seems to correlate with smaller deviation 