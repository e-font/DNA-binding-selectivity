import numpy as nm
import numpy.random as rnd
import model1
import correlations
import matplotlib.pyplot as plt
A = 1
bases = ['A', 'T', 'C', 'G']
pairs = nm.array([[p + q for q in bases] for p in bases])
energy = rnd.normal(0.5, 0.05, (16, 16))
energy = rnd.uniform(0, 0.5, (16,16))
nm.fill_diagonal(energy, 0)
energy = energy + energy.T
#print energy


def pair_energy(pair1, pair2):
    def pair_index(pair):
        return bases.index(pair[0]) * 4 + bases.index(pair[1])
    return energy[pair_index(pair1), pair_index(pair2)]

def partition_function(x, tar_seq, gen_seq):
    z = 0
    for i in range(0, len(gen_seq) - len(tar_seq)):
        e = 0.
        for j in range(len(tar_seq)/2 + 1):
            p1 = tar_seq[j] + tar_seq[j + 1]
            p2 = gen_seq[i + j] + gen_seq[i + j + 1]
            e += pair_energy(p1, p2)
        z += nm.exp(-1 * x * e)
    return z

def q_perfect_match(x, tar_seq, gen_seq):
    q = 0
    for i in range(0, len(gen_seq) - len(tar_seq)):
        if all(tar_seq[j] == gen_seq[j + i] for j in range(len(tar_seq))):
            q += 1
    return q

def factorized_z(x, tar_seq, p_matrix):
    z = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        z_b = 0.
        for j, row in enumerate(pairs):
            for k, p in enumerate(row):
                # print nm.log(p_matrix[j, k])
                if pair == p:
                    z_b += p_matrix[j, k]
                else:
                    z_b += p_matrix[j, k] * nm.exp(-1. * x * pair_energy(p, pair))
        print z_b
        z = z * z_b
    return z

def factorized_q(tar_seq, p_matrix):
    q = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        z_b = 0.
        for j, row in enumerate(pairs):
            for k, p in enumerate(row):
                if pair == p:
                    z_b = p_matrix[j, k]
                    break
        print z_b
        q = q * z_b
    return q

if A == 1:
    p_bases = [0.25, 0.25, 0.25, 0.25]
    P = rnd.normal(1, 0.1, (4, 4))
    P = P / nm.sum(P)
    l = 4
    N_g = 100000

    x = nm.linspace(0, 10)
    t = rnd.choice(bases, size=l, p=map(lambda x: float(x), p_bases))
    g = ''
    for i in range(N_g / 2):
        g += correlations.random_pair(P)

    Y = lambda x: model1.selectivity(partition_function(x, t, g), q_perfect_match(x,t,g))
    y = lambda x: model1.selectivity(factorized_z(x, t[1:0:1], P) + factorized_z(x, t, P), factorized_q(t[1:0:1], P) + factorized_q(t, P))
    #factorized_q(t, P)
    #factorized_z(2, t, P)
    #print partition_function(10,t,g), q_perfect_match(10, t, g)
    print factorized_z(x, t[1:0:1], P), factorized_z(x, t, P), factorized_q(t[1:0:1], P), factorized_q(t, P)
    plt.plot(x, Y(x), x, y(x))
    plt.show()