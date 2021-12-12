import numpy as np
import matplotlib.pyplot as plt

np.random.seed(100000)
data = np.random.uniform(0.0,1.0,100000)

plt.hist(data)
plt.show()
def gaussa(r1,r2):
    return np.sqrt(-2*np.log(r1))*np.cos(2*np.pi*r2)
def gaussb(r1,r2):
    return np.sqrt(-2*np.log(r2))*np.sin(2*np.pi*r1)
gauss = np.zeros(100000)
for i,(a,b) in enumerate(zip(data[::2],data[1::2])):
    idx = 2*i
    gauss[idx] = gaussa(a,b)
    gauss[idx+1] = gaussb(a,b)
plt.hist(gauss,bins="auto")
plt.show()