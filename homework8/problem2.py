import scipy
import scipy.special
import numpy as np

N = 50
M = N//6

def direct(m):
    def fac(k):
        return scipy.special.factorial(k,exact=False)
    #return choose(N,2) / (choose(m,2)*choose(N-m, 2))
    return fac(N)/(fac(m)*fac(N-m))
    #computing this directly is so difficult because these numbers are so large
    #thus, they cost precision to represent and the error multiplies itself
    #very quickly

def coincidences_approx(m):
    #returns estimate of the number of ways exactly m particles out of N can be excited
    iterations = 100000
    def gen_state(): # represents the energies of all the particles
        bits = np.random.choice(range(N),size=m,replace=False)
        asnum = 0
        for x in bits:
            asnum = asnum | (1 << x)
        return asnum # each is one of many many possible states 

    coin = dict()
    for _ in range(iterations):
        s = gen_state()
        if s in coin:
            coin[s] += 1 
        else:
            coin[s] = 1 
    num_pairs_of_states = lambda num_states: num_states*(num_states-1)/2
    ss = 0
    for x in coin.values():
        ss += num_pairs_of_states(x)
    return num_pairs_of_states(iterations)/ss
    #I understand it now!



coincidences = np.zeros(M)
coincidences2 = np.zeros(M)
for m in range(2,M):
    coincidences[m] = direct(m)
    #the second method, outlined in the hint
    coincidences2[m] = coincidences_approx(m)

import matplotlib.pyplot as plt

plt.plot(np.array(range(M)), coincidences, label="direct formula")
plt.plot(np.array(range(M)), coincidences2, label="indirect formula")
plt.legend()
plt.title("N="+str(N))
plt.xlabel("m")
plt.ylabel("coincidences")
#plt.yscale("log")
plt.show() #it seems that at ~ 5.5 the formula becomes unstable, 
#as that is when it diverges from the coincedences approach