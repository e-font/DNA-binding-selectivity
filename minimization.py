import numpy as nm
import scipy as spi
import matplotlib.pyplot as plt

a = 0.5
N = 1.

def h(x):
    return a*x - nm.log(nm.c)

def f(x):
    return a*x + x*nm.log(x) + (N - x) * nm.log(N - x)

def g(x):
    return f(x) + 0.5 * nm.log(2. * nm.pi * x) + 0.5 * nm.log(2. * nm.pi * (N - x))

x = nm.linspace(0., 1., 100)

plt.plot(x, g(x))
plt.show()