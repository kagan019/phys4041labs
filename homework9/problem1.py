import numpy as np
import matplotlib.pyplot as plt

# Peter,
# What book did the chapter Professor Vinals linked for us come from?

#task:
#use monte carlo to find equilibrium/canonical distribution of 2D ising model

#assume arbitrary units so kBT = 1 and J = J' [unitless], numerically speaking
Js = [0.4,0.5] 
steps = 10000
N= 100 #NXN lattice sites
dE = 0.01

def initialize_lattice():
    lat = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            lat[i][j] = int(np.random.choice([-1,1]))
    return lat

def flip_random_dipole(lattice):
    lattice2 = [[x for x in y] for y in lattice]
    i = np.random.randint(0,N)
    j = np.random.randint(0,N)
    lattice2[i][j] = [1, None, -1][int(lattice2[i][j] + 1)]
    return lattice2

def energy_of_configuration(lattice): #unitless (numerically E/(kB*T))
    #adjacent dipoles 'want' to be like one another
    #periodic boundaries; torus topology
    s = 0
    for i in range(N-1): # vertical pairs
        for j in range(N):
            s += lattice[i][j] * lattice[i+1][j]
    #periodic boundary conditions
    for j in range(N):
        s += lattice[N-1][j] * lattice[0][j]
    for i in range(N):
        for j in range(N-1): # horizontal pairs
            s += lattice[i][j] * lattice[i][j+1]
    for i in range(N):
        s += lattice[i][N-1] * lattice[i][0]
    return -J * s

def spinsum(lattice):
    s = 0
    for i in range(N):
        for j in range(N):
            s += lattice[i][j]
    return s

def metropolis(S1,S2): #likelihood of S1 -> S2
    energyS1 = energy_of_configuration(S1)
    energyS2 = energy_of_configuration(S2)
    deltaE = energyS2 - energyS1
    return 1 if energyS2 <= energyS1 else np.exp(-deltaE)

figures = 0
def ising(J):
    global figures
    energy_bins = np.linspace(-J*N*N,J*N*N,int(2*J*N*N/dE))
    def bin_config(lattice): 
        egy = energy_of_configuration(lattice)
        idx_of_egy = int((egy + J*N*N) / dE)
        if idx_of_egy < 0 or idx_of_egy > len(energy_bins):
            return None
        return idx_of_egy
    occurence_of_energy = np.zeros(len(energy_bins))
    current_lattice = initialize_lattice()  
    
    for step in range(steps):
        consider = flip_random_dipole(current_lattice)
        die = np.random.uniform()
        if die < metropolis(current_lattice,consider):
            current_lattice = consider
        idx = bin_config(current_lattice)
        assert(idx is not None)
        occurence_of_energy[idx] += 1

                                                     #  J     0.4       0.5  
    print(energy_of_configuration(current_lattice)/N**2) #  ~-0.2     ~-0.3
    print(1/N**2*spinsum(current_lattice))               #   ~0.0      ~0.0
    plt.figure(figures)
    figures += 1
    largest_num_occurences = max(occurence_of_energy) 
    relative_frequency = [x/largest_num_occurences for x in occurence_of_energy]
    plt.plot(energy_bins,relative_frequency)
    plt.title("J=" + str(J))
    plt.xlabel("energy (kBT)")
    plt.ylabel("relative frequency") #should appear like canonical distribution of ising model
    

for J in Js:
    ising(J)
plt.show()