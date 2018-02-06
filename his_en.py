"""
Creates an histogram for the Santa Lucia distribution of energies from a matrix.
"""

import numpy as nm
import matplotlib.pyplot as plt
a = nm.loadtxt("santa_lucia04.txt") # The matrix
a = nm.rot90(a)
b = nm.array(nm.diag(a))
nm.fill_diagonal(a, 0.)
plt.hist([i for i in a.flatten() if i != 0.], label='x')
plt.hist(b, label='y')
plt.show()
