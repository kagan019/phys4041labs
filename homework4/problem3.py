import numpy as np
import matplotlib.pyplot as plt 
import time
np.random.seed(int(time.time()))

maxt = 10000
dt = 0.05
lsp = np.linspace(0,maxt,num=int(maxt/dt))

magb1 = (0.1,0.1,0.1,0.38)
magb2 = (0.1,0.1,0.1,0.41)
def R0(magb): 
    m,a,g,b = magb
    return a/(m+a)*b/(m+g)
def rand_init_SEIR(N):
    left = N
    ret = []
    for _ in range(3):
        chosen = 0 if left == 0 else np.random.randint(0,left)
        ret += [chosen]
        left -= chosen
    ret += [N-sum(ret)]
    return ret

# dS/dt = s(t)
# dE/dt = e(t)
# dI/dt = i(t)
# dR/dt = r(t)
def s(SEIR, magb):
    N = sum(SEIR)
    S,E,I,R = SEIR
    m,a,g,b = magb
    return m*N-m*S-b*I*S/N
def e(SEIR, magb):
    N = sum(SEIR)
    S,E,I,R = SEIR
    m,a,g,b = magb
    return b*I*S/N-(m+a)*E
def i(SEIR, magb):
    N = sum(SEIR)
    S,E,I,R = SEIR
    m,a,g,b = magb
    return a*E-(g+m)*I
def r(SEIR, magb):
    N = sum(SEIR)
    S,E,I,R = SEIR
    m,a,g,b = magb
    return g*I-m*R
def next_SEIR(SEIR, magb):
    # Runge-Kutta order 2
    K1 = [dt*(f(SEIR,magb)) for f in [s,e,i,r]]
    SEIRmid = [y+k1/2 for y,k1 in zip(SEIR,K1)]
    K2 = [dt*(f(SEIRmid,magb)) for f in [s,e,i,r]]
    return [y+k2 for y,k2 in zip(SEIR, K2)]

#run simulation
sims = {}
for N in range(200,1001,200):
    # 3 initial conditions
    print(N)
    init = rand_init_SEIR(N)
    SEIRs_by_coeff = [None, None]
    for j,magb in enumerate((magb1,magb2)):
        SEIR_evol = []
        for jj,ti in enumerate(lsp):
            SEIR_evol += [
                init if jj == 0 
                else next_SEIR(SEIR_evol[jj-1], magb)
            ]
        SEIRs_by_coeff[j] = SEIR_evol
    sims[N] = SEIRs_by_coeff

#format plots
_,ax =plt.subplots(4,2, sharex=True)
yaxis_names = ["susceptible","exposed","infected","recovered"]
for bc in range(2):
    for dim in range(4):
        myind = (dim,bc)
        ax[myind].set_ylabel("number {}".format(yaxis_names[dim]))
plt.xlabel("log of time")
plt.xscale("log")
#plot points
for k in sims:
    s = sims[k]
    for bc, bycoeff in enumerate(s):
        to_plot_individually = np.transpose(bycoeff) #a list of 4-tuples to a 4-tuple of lists
        for sim_dimension in range(4): #S, E, I,or R 
            myind = sim_dimension,bc
            ax[myind].set_title("magb coeffs {}".format(bc+1))
            ax[myind].plot(lsp,to_plot_individually[sim_dimension], label="N = {}".format(k))
            ax[myind].legend(loc="upper right")

print("R0 parameter {} and {} for magb coeffs 1 and 2, respectively".format(R0(magb1), R0(magb2)))
#R0 parameter 0.95 and 1.025 for magb coeffs 1 and 2, respectively

plt.show()
# indeed, all curves behave asymptotically.
#lets pray our curve is a lucky one