import numpy as nm
import scipy.special as spe
import matplotlib.pyplot as plt

def monomial(x, c, n, i):
    if i % 2 == 0:
        j = i / 2.
        return (-1) ** j * x ** (n - 2 * j) * c ** j * spe.binom(n - j, j)
    else:
        return 0

def polynomial(x, c, n):
    return [monomial(1, c, n, j) for j in range(0, n + 1)]

a = 1
b = 1
c = -1. * (a + b)

x = nm.arange(1, 40)
for xe, ye in zip(x, [c - nm.roots(polynomial(c, a * b, j)) for j in x]):
    plt.scatter([xe] * len(ye), ye)
    
plt.show()