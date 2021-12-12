#not totally functional, but the ideas are there
#the issue is the stability of the van Beijeren function.
#check out the implementation below.

import numpy as np
from scipy import stats

T = 1 #units of ep/k_B 
m = 1 # nuanced; simultaneously, the unit of mass, 
      #and, the unitless quantity representing how many such mass units constitute each particle
beta = 1/(1*T) #1/k_B*T #the unit of inverse energy
#beta = epd/ep = 1/ep for epd=1, so ep=k_BT
Ld = 10 #dimensionless side length of bounding box
N = 100 #num of particles
sigma = 1 #arbitrary length units
tau = np.sqrt(m*sigma**2*beta)  #the chosen unit of time
#differentials
dt = tau/1000
assert(m*(sigma/tau)**2 - 1/beta < 0.00001) #energy units (ep/(k_B*T) = 1)

def mag(vec):
    return np.sqrt(vec[0]**2+vec[1]**2)

class Particle:
    def __init__(self,r,v):
        self.r = r #dimensionless count of sigma
        self.v = v #dimensionless count of sigma/tau

def thermometer(particle):
    def van_Beijeren(v): #unitless v
        #v is the normal component of the particle velocity
        assert(v != 0)
        #ep/(k_B*T) = 1
        return np.sqrt(-2*np.log(1-np.exp(-v**2/2))) #dimensionless vel

    if particle.r[0] >= Ld and particle.v[0] > 0:
        particle.v[0] = -van_Beijeren(particle.v[0])
        particle.r[0] = Ld
    elif particle.r[0] <= 0 and particle.v[0] < 0:
        particle.v[0] = van_Beijeren(particle.v[0])
        particle.r[0] = 0
    if particle.r[1] >= Ld and particle.v[1] > 0:
        particle.v[1] = -van_Beijeren(particle.v[1])
        particle.r[1] = Ld
    elif particle.r[1] <= 0 and particle.v[1] < 0:
        particle.v[1] = van_Beijeren(particle.v[1])
        particle.r[1] = 0
    

#potential function
def V(r): # r unitless
    if r == 0:
        exit("too close!")
    # V(r) = 4/beta*((sigma/r)**12-(sigma/r)**6), units of 1/beta, r units of sigma
    # V, converted to unitless.
    return 4*((1/r)**12-(1/r)**6)

def dVdr(r): # r unitless
    if r == 0:
        exit("too close!")
    # dV'/dr' (unitless).
    return 4*((-12)/r**13-(-6)/r**7)

def forces(particles): 
    #return the total force on each particle as a vector
    f = np.zeros((len(particles),2))
    for i in range(len(particles)-1):
        for j in range(i+1,len(particles)):
            pi = particles[i]
            pj = particles[j]
            a = pi.r-pj.r #unitless, because r is unitless
            duijdr = dVdr(mag(a)) #unitless
            df = a/mag(a)*duijdr #still unitless, but with direction
            f[i] -= df
            f[j] += df
    return f #unitless force

def velverlet(particles):
    for p in particles:
        thermometer(p)

    fs = forces(particles) #unitless
    for i,p in enumerate(particles):
        #units of sigma times the unitless expression
        p.r = p.r + (1/tau)*dt*p.v+(m/tau**2)*dt**2/(2*m)*fs[i]
        #(here, mass cancels out because the mass of the particle 
        # equals the unit mass)

    fs2 = forces(particles)
    for i,p in enumerate(particles):
        #units of sigma/tau times the unitless expression
        p.v = p.v + (m/tau)*dt/(2*m)*(fs2[i]+fs[i])

def instantiate_fluid():
    #according to bolztmann distribution of velocities in a lattice
    #remember this is a 2 dimensional simulation so there are 2, not 3,
    #degrees of freedom
    ptcls = np.full(N,None)
    for i in enumerate(ptcls):
        ptcls[i] = Particle(np.zeros(2),np.zeros(2))
    L = Ld*sigma #units of sigma
    numside = int(np.sqrt(N))
    xivl = L/numside
    yivl = L/numside
    xpos = np.linspace(0,L-xivl,numside) +xivl/2
    ypos = np.linspace(0,L-yivl,numside) +yivl/2
    for i,x in enumerate(xpos):
        for j,y in enumerate(ypos):
            idx = i*numside+j
            # particle positions are unitless counts of sigma
            ptcls[idx].r = np.array([x,y]) / sigma 
    #distribute velocities
    #from Wikipedia:
    # Mathematically, the Maxwell–Boltzmann distribution is the chi distribution with three degrees of freedom (the components of the velocity vector in Euclidean space), with a scale parameter measuring speeds in units proportional to the square root of {\displaystyle T/m}{\displaystyle T/m} (the ratio of temperature and particle mass).[2]
    # from Scipy docs:
    # A special case of a chi distribution, with df=3, loc=0.0, and given scale = a, where a is the parameter used in the Mathworld description [1].
    a = np.sqrt(1/(beta*m)) #units of sigma/tau
    mxwl = stats.chi(scale=a, loc=0.0, df=2) # 2 dof for 2d space!
    speeds = mxwl.rvs(size=N) #units of sigma/tau, due to the units of the scale parameter
    angles = 2*np.pi*np.random.rand(N) #unitless (or radians)
    for i,p in enumerate(ptcls):
        p.v = np.array([speeds[i]*np.cos(angles[i]), 
                    speeds[i]*np.sin(angles[i])]) / (sigma/tau) #unitless count of sigma/tau

    return ptcls

def ptcl_energy(particles):
    vsq = 0
    for p in particles:
        vsq += mag(p.v)*mag(p.v)
    totV = 0
    for i,pi in enumerate(particles[:-1]):
        for pj in particles[i+1:]:
            totV += V(mag(pi.r-pj.r))
    return 1/2*m*vsq + totV #unitless

# macroscopic metrics
def avg_interptcl_sep(particles):
    ravg = 0
    for i,pi in enumerate(particles[::-1]):
        for j,pj in enumerate(particles[i::]):
            ravg += mag(particles[i].r-particles[j].r)
    
    return ravg/int(N*(N-1)/2) #unitless

def avg_velsq(particles):
    vsq = 0
    for p in particles:
        vsq += mag(p.v)*mag(p.v)
    return vsq/N #unitless

# simulator
steps = 1000
print(steps)
times = np.zeros(steps)
Ut = np.zeros(steps)
Ravgt = np.zeros(steps)
vsqavg = np.zeros(steps)
particles = instantiate_fluid()
assert(len(particles) == N)

t = 0
for i in range(steps):
    if i % 100 == 0:
        print(i) 
    times[i] = t
    Ut[i] = ptcl_energy(particles)
    Ravgt[i] = avg_interptcl_sep(particles)
    vsqavg[i] = avg_velsq(particles)
    t += dt
    velverlet(particles)

import matplotlib.pyplot as plt

plt.plot(times/tau,Ut,label="internal energy")
plt.plot(times/tau,Ravgt,label="avg interparticle separation")
plt.plot(times/tau,vsqavg,label="avg velocity squared")
plt.legend()
plt.xlabel("dimensionless time")
plt.ylabel("dimensionless E, r, or v^2")
plt.show()


