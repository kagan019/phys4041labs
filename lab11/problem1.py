import numpy as np
import matplotlib.pyplot as plt

A = 0.5
B = 0.2
X = np.linspace(0,10)
plt.plot(X, np.sqrt(B/(2*np.pi*A))*np.exp((-B/(2*A))*(X-A/B)**2))

#state vector
n = 100
t = 0

def deltat():
    return -1/(A+B*n)*np.log(np.random.uniform())

bins = np.zeros(1000)
while True:
    r = np.random.uniform()
    if r < A/(A+B*n):
        n += 1
    elif r > A/(A+B*n):
        if n != 0:
            n -= 1
    t += deltat()
    bins[n] += 1

tot = sum(bins)

dp = 0.01

relbins = np.zeros(0,10,int(10/dp))
for i,count in enumerate(bins):
    relprb = i/tot
    idx = int(relprb/dp)
    relbins[idx] = relprb

horz = np.linspace(0,10,10)
np.plot(horz,relbins)
np.show()

