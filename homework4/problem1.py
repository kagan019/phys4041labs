#p (mmHg)	0.8	0.4	0.2	0.1	0.05
#x	        740	487	475	485	489

# where x(p=0) is the value I ultimately wish to compute.
# I begin by writing
#       x=c0+c2*p**2+c3*p**3+c6*p**6,  ci != 0
# so the repeated Richardson extrapolation applies as if p was
# the grid spacing of a numerical computation.

# Call the list of exponents 
ss = [None, 2, 3, 6]
# The value p = 0.8/2**c is known for each c in [0, 1, 2, 3, 4]
def x(c): #x(p/2**c)
    return [740,487,475,485,489][c]

A = [[x(m)] for m in range(len(ss))]
for m in range(len(ss)):
    for k in range(1,m+1):
        A[m] += [A[m][k-1] + (A[m][k-1]-A[m-1][k-1])/(2**(ss[k])-1)]
print(A[len(ss)-1][len(ss)-1])
print(str(A))