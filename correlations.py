import numpy as nm
import model1
import matplotlib.pyplot as plt

pairs  = 5.

n_2 = 5.
n_1 = 0.
n_0 = pairs - n_2 - n_1

d = -0.8

def z_0(x):
    return (0.25 * (nm.exp(x) + 3)) ** (pairs * 2.)

def q_0(x):
    return (0.25 * (nm.exp(x))) ** (pairs * 2)

def z(x):
    return z_0(x) * (1. + (d/15.) * (16. * nm.exp(2*x) / ((nm.exp(x) + 3) ** 2) - 1)) ** n_2 \
           * (1. + (d/15.) * (16. * nm.exp(x) / ((nm.exp(x) + 3) ** 2) - 1)) ** n_1 \
            * (1. + (d/15.) * (16. / ((nm.exp(x) + 3) ** 2) - 1)) ** n_0

def q(x):
    return q_0(x) * (1. + d) ** n_2 * (1. - (d/15.)) ** (n_1 + n_0)

x = nm.linspace(0, 0.2)
y_0 = nm.vectorize(lambda i: model1.selectivity(z_0(i), q_0(i)))
y = nm.vectorize(lambda i: model1.selectivity(z(i), q(i)))

plt.plot(x, z_0(x),x,  z(x))
plt.show()
