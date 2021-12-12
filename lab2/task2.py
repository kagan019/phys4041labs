import math

x_values = [0.10,1.00,10.0]
selected_iters = [3,5,8]

def J0(x):
    return math.sin(x)/x

def J1(x):
    return (math.sin(x) - x*math.cos(x))/x**2

budp = []
def besselup(x : float):
    global budp
    budp = [0] * (max(selected_iters)+1)
    budp[0] = J0(x)
    budp[1] = J1(x)
    i = 0
    while i <= max(selected_iters)-2:
        coeff = (2*(i+1)+1)/x
        budp[i+2] = coeff*budp[i+1]- budp[i]
        i += 1
bddp = []
def besseldown(x : float):
    global bddp
    bddp = [0] * (max(selected_iters)+1)
    bddp[0] = J0(x)
    bddp[1] = J1(x)
    i = max(selected_iters)*10
    n, c, v = 0,0.5,0.5
    while i >= 2:
        coeff = (2*(i-1)+1)/x
        n = coeff*c - v
        if i < len(bddp):
            bddp[i-2] = n
            bddp[i-1] = c
            bddp[i] = v
        v = c
        c = n
        i -= 1
    scale = J0(x)/bddp[0]
    #print("scale: {}".format(scale))
    for i,_ in enumerate(bddp):
        bddp[i] *= scale

def test_bessel(besselfun, up=True):
    print("Bessel J_l(x) at various iterations l by input value x")
    print("x\t{}\t\t{}\t\t{}".format(*selected_iters))
    for x in x_values:
        besselfun(x)
        print("{}\t".format(x),end="")
        for i in selected_iters:
            print("{x:.6e}".format(x=budp[i] if up else bddp[i]), end = "\t")
        print() 
print("UPWARDS")
test_bessel(besselup, up=True)
print("DOWNWARDS")
test_bessel(besseldown, up=False)