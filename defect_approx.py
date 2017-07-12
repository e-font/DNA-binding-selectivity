import numpy as nm
import scipy.optimize as op
import matplotlib.pyplot as plt

def f0(x, N, i):
    return nm.log(x) - nm.log(N - x) - i


def f1(x, N, i):
    return f0(x, N, i) + 0.5 * ((x) ** -1. - (N - x) ** -1.)


def a0(N, i):
    return N * nm.exp(i) / (1. + nm.exp(i))


def a1(N, i):
    return a0(N, i) * (2. * N * nm.exp(i) - 2.) / (2. * N * nm.exp(i) - nm.exp(i) - 1.)

N = 10.
steps = 1000
c = nm.linspace(-1., 1., num=50)
# x = nm.linspace(steps ** -1., N, num=steps, endpoint=False)
x = nm.linspace(0.5, N - 0.5, num=1000)

#print f1(N - 0.001, 10., 0.)
y0 = nm.vectorize(lambda i: op.brenth(f0, 0.5,  N - .5, args=(N, i)))(c)
y1 = nm.vectorize(lambda i: op.brenth(f1, 0.5,  N - .5, args=(N, i)))(c)
w0 = nm.vectorize(lambda i: a0(N, i))(c)
w1 = nm.vectorize(lambda i: a1(N, i))(c)

# op.brenth(f0, 0, N, args=(N, 3))
#plt.plot(x, f0(x, N, -1.), x, f1(x, N, -1.))
plt.plot(c, y1 - w1)
# plt.plot(x, f1(x, N, 0.))

plt.show()