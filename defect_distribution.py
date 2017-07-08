import numpy as nm
from numpy import random as rnd
import matplotlib.pyplot as plt
import model1

def defect_count(src_seq, tar_seq):
    hams = nm.array(model1.hamiltonians(src_seq, tar_seq))
    lst = []
    for i in range(src_seq.size):
        lst.append((hams[:,1] == i).sum())
    return lst


def defect_count_base(base, src_seq, tar_seq):
    s_ary = src_seq
    s_ary[s_ary != base] = 'N'
    hams = nm.array(model1.hamiltonians(s_ary, tar_seq))
    hams[:, 1] -= (s_ary != base).sum()
    lst = []
    for i in range(src_seq.size):
        lst.append((hams[:, 1] == i).sum())
    return lst


if __name__ == '__main__':
    # BASES = nm.array([['A', 'T', 'C', 'G'], [0.2922, 0.2928, 0.2074, 0.2076]])
    BASES = nm.array([['A', 'T', 'C', 'G'], [0.25, 0.25, 0.25, 0.25]])
    beta = 10
    n = 30
    N = 100000
    E_MATCH = -1.
    E_NO_MATCH = 0.

    source_sequence = rnd.choice(BASES[0], size=n, p=map(lambda x: float(x), BASES[1]))
    target_sequence = rnd.choice(BASES[0], size=N, p=map(lambda x: float(x), BASES[1]))

    occurrences = map(lambda b: sum(b == j for j in source_sequence), BASES[0])
    print occurrences[0]

    defects = nm.array(defect_count_base('A', source_sequence, target_sequence))
    x = nm.arange(n)
    y = defects * nm.exp(-1. * beta * (E_MATCH * (n - x) + E_NO_MATCH * x))
    q = y[0]
    z = y.sum()
    avg = n * nm.average(x * y / z)

    # print y[round(avg)], z
    print avg, defects

    # plt.plot(x, y)
    # plt.show()