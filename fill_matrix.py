import numpy as nm
bases = ['A', 'C', 'G', 'T']
compl_bases = ['T', 'G', 'C', 'A']
pairs = nm.array([[p + q for q in bases] for p in bases])
energies = nm.zeros((16,16))
for i in range(16):
    for j in range(16):
        if j/4 == (3 - i/4):
            energies[i][j] = 1
energies1 = nm.zeros(energies.shape)
energies = nm.loadtxt("santa_lucia04.txt", delimiter='\t')
def index_to_base(i):
    return bases[i]

def base_to_index(X):
    return bases.index(X)

def index_to_pair(i):
    return pairs[i / 4][i % 4]

def pair_to_index(XY):
    for i, p in enumerate(pairs.flatten()):
        if p == XY:
            return i
    raise Exception

def what_value(I, J):
    i = int(I)
    j = int(J)
    if not energies[i][j] == 0.:
        print energies[i][j]
        return energies[i][j]
    else:
        sum1 = sum(what_value(k * 4 + (i % 4), (k + 1) * -4. + (j % 4))
                   for k in range(4))
        sum2 = sum(what_value((i / 4) * 4 + k, (j / 4) * 4 - (k + 1))
                   for k in range(4))
        return (sum1 + sum2) / 8.

def fill():
    for i in range(16):
        for j in range(16):
            energies1[i][j] = what_value(i, j)
fill()
print energies1
#nm.savetxt('energies1.txt', energies1)
