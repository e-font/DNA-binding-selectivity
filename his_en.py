import numpy as nm
import matplotlib.pyplot as plt
a = nm.loadtxt("santa_lucia04.txt")
a = nm.rot90(a)
b = nm.diag(a)
plt.hist([i for i in (a - b).flatten() if i != 0.], label='x')
plt.hist([i for i in b.flatten() if i != 0.], label='y')
plt.show()