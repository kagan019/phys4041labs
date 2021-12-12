import matplotlib.pyplot as pyplot
import math

def test(x):
  return math.sin(x)/(x*x)-math.cos(x)/x

def bessel(x,n):
  i = 0
  f = 1
  s = 0
  while i <= n:
    v = pow(x,i)/f
    v /= x
    if i % 2 == 1:
      v /= x
    if (i+1)% 4 < 2: 
      v = -v
    s += v 
    i += 1
    f *= i
  return s
x_ = [0.02]
b = []
while x_[-1] <= 10.0:
  b += [bessel(x_[-1],80)]
  x_ += [x_[-1] + 0.02]
pyplot.plot(x_[:-1],b)
pyplot.show()

