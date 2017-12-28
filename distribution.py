import numpy as nm
import numpy.random as rnd
import matplotlib.pyplot as plt

bases = [0, 1, 2, 3]
p_distr = nm.full((4,), 1 / 4.)

energies = rnd.normal(0.5, 0.1, (16, 16))
energies = energies + energies.T
#energies = nm.full((16, 16), 1.)
nm.fill_diagonal(energies, 0)
energies[2 * 4 + 1][2 * 4 + 1] = -.1  # Corresponds to (2,1) = GC pair
#energies = nm.rot90(nm.loadtxt("energies1.txt"))

def random_seq(l):
    return rnd.choice(bases, size=l, p=p_distr)


def part_func(t_seq):
    z = 1.
    for i in range(0, len(t_seq), 2):
        pair = t_seq[i] * 4 + t_seq[i + 1]
        e = 0.
        for j in range(16):
            e += nm.exp(-10. * energies[pair][j])
        z *= e
    return z


def q_match(t_seq):
    z = 1.
    for i in range(0, len(t_seq), 2):
        pair = t_seq[i] * 4 + t_seq[i + 1]
        e = nm.exp(-10. * energies[pair][pair])
        z *= e
    return z


def selectivity(t_seq, step_size=2):
    z = part_func(t_seq)
    #q = q_match(t_seq)
    if step_size == 1:
        z += part_func(t_seq[1:-1])
        #q += q_match(t_seq[1:-1])
    return z


def entropy(t_seq):
    s = 0.
    l = len(t_seq)
    print list(t_seq).count(1)
    occurrences = [float(list(t_seq).count(i)) / l for i in bases]
    for o in occurrences:
        if not o == 0:
            s += o * nm.log(o)
    return -1 * s


s = []
t = nm.array([0, 1,2,2, 2, 2, 2,3,3,3,3,3, 3, 3, 3, 3])
t1 = (t + 1) % 4
t2 = (t + 2) % 4
t3 = (t + 3) % 4
for i in range(10 ** 4):
    s.append(nm.log(selectivity(rnd.permutation(t), step_size=2)))
    s.append(nm.log(selectivity(rnd.permutation(t1), step_size=2)))
    s.append(nm.log(selectivity(rnd.permutation(t2), step_size=2)))
    s.append(nm.log(selectivity(rnd.permutation(t3), step_size=2)))

# nm.save("10^6_t_20_RAW", s)
# print sorted(s, reverse=True)
#h = nm.histogram(s, bins=100)
#print h[1]
# nm.savetxt("10^5_t_20_perm_0001_anomaly.csv", h[0], delimiter=';')
plt.ylabel('Frequency')
plt.xlabel('ln selectivity')
plt.title('Permutated array, step = 1. (Real energies)')
plt.hist(s, bins=100)
plt.show()

