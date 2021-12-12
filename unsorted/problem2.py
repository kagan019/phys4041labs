import numpy as np
import matplotlib.pyplot as plt

def a0(q):
    return -1/2*q**2+7/128*q**4-29/2304*q**6
print(a0(1))
def b1(q):
    return 1-q-1/8*q**2+1/64*q**3-1/1536*q**4
print(b1(1))

for init in [1,2]:
    q=1
    p=a0(q)
    sysinit = [1.0, 0.0]
    if init == 2:
        p=b1(q)
        sysinit = [0.0, 1.0]

    N = 10000
    tspace = np.linspace(0,30,N)
    dt = tspace[1]-tspace[0]

    # x = x(t)
    # dx/dt = z(t)
    # dz/dt = -(p-2*q*cos(2*t))*x
    def f(t,sys):
        x,z = sys
        dzdt = -(p-2*q*np.cos(2*t))*x
        return np.array([z,dzdt]) 

    def rko4():
        sys = np.array([sysinit[0],sysinit[1]])
        for t in tspace:
            yield (t,sys[0])
            k1 = dt*f(t,sys)
            k2=dt*f(t+dt/2,sys+k1/2)
            k3=dt*f(t+dt/2,sys+k2/2)
            k4=dt*f(t+dt, sys+k3)
            sys += 1/6*(k1+2*k2+2*k3+k4)

    X, Y = [],[]
    for (x,y) in rko4():
        X += [x]
        Y += [y]
    plt.plot(X,Y, label="p=a0(q)" if init == 1 else "p=b1(q)")
plt.legend()
plt.xlabel("t")
plt.ylabel("x")
plt.show()
#discussion:
# the a0 plot average amplitude is described roughly by apparent exponential growth, so it is continually resonating
# the b1 plot is fairly balanced between damping and resonance. 