"""
Main script for model3. Contains all functions necessary to calculate the computed and the actual selectivity
Match_assured: set true is equivalent to adding a matching string to a position in the genome sequence. Default is false.
Step_size: the number of steps the targeting sequence 'slides' on the genome one. Default is 2 for simplicity, real world is 1.

"""

import numpy as nm
import numpy.random as rnd
import model1
import correlations
import matplotlib.pyplot as plt

bases = ['A', 'C', 'G', 'T']
pairs = nm.array([[p + q for q in bases] for p in bases])  # matrix of base pairs

def compl_seq(seq):
    c = ''
    for b in seq:
        i = int(-1. * bases.index(b)) - 1
        c += bases[i]
    return c

def pair_index(pair):  # Returns the index of the base duplex in the energy matrix
    return bases.index(pair[0]) * 4 + bases.index(pair[1])

def random_pair(P):  # Generates a random duplex given a probability distribution for all the duplexes P
    i = nm.random.choice(range(4), 1, p=nm.sum(P, axis=1))[0]
    j = nm.random.choice(range(4), 1, p=P[i] / nm.sum(P[i]))[0]
    return pairs[i, j]

def random_seq(size, p_bases):  # Generates a random targeting sequence given a probability distribution (not by duplexes)
    return rnd.choice(bases, size, p=map(lambda x: float(x), p_bases))

def random_seq_gen(size_gen, P):  # Generates a random genome sequence given a length and a probability distribution P (by duplexes)
    p = rnd.choice(pairs.flatten(), size_gen / 2, p=map(lambda x: float(x), P.flatten()))
    return ''.join(p)

def pair_energy(pair1, pair2):  # Returns binding energy given two base pairs
    return energy[pair_index(pair1), pair_index(pair2)]

def partition_function(beta, tar_seq, gen_seq, step_size=2, match_assured=False):  # Partition function, computed with the actual sequences
    z = 0.
    i = 0  # Counts the position within the genome sequence
    while i <= len(gen_seq) - len(tar_seq):
        e = 0. # total energy for a given position of tar seq along gen seq
        for j in range(0, len(tar_seq), 2):
            p1 = tar_seq[j] + tar_seq[j + 1]
            p2 = gen_seq[i + j] + gen_seq[i + j + 1]
            e += pair_energy(p1, p2)
        z += nm.exp(-1. * beta * e)
        i += step_size  # Steps by step_size
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        z += nm.exp(-1. * beta * e)
    return z

def q_perfect_match(beta, tar_seq, gen_seq, step_size=2., match_assured=False):  # Actual boltzmann probability of perfect matches
    q = 0.
    i = 0
    while i <= len(gen_seq) - len(tar_seq):
        if all(tar_seq[j] == gen_seq[j + i] for j in range(len(tar_seq))):  # Extra condition compared to before: all bases must match
            e = 0.
            for j in range(0, len(tar_seq), 2):
                p1 = tar_seq[j] + tar_seq[j + 1]
                p2 = gen_seq[i + j] + gen_seq[i + j + 1]
                e += pair_energy(p1, compl_seq(p2))
            q += nm.exp(-1 * beta * e)
        i += step_size
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        q += nm.exp(-1. * beta * e)
    return q

def factorized_z(beta, tar_seq, len_gen_seq, p_matrix, step_size=2, match_assured=False):  # Mathematical model of Z, factorized. Mean field value
    z = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        z_b = 0.  # Z for the current base pair
        for j, row in enumerate(pairs): # running over all 16 possible base pairs
            for k, p in enumerate(row):
                z_b += p_matrix[j, k] * nm.exp(-1. * beta * pair_energy(pair, p))
        z = z * z_b
    z = z * (len_gen_seq - len(tar_seq) + 1) / 2.  # Rescales Z so that it can be compared with the computed one
    if step_size == 1:
        z += factorized_z(beta, tar_seq, len_gen_seq, p_matrix, step_size=2, match_assured=False)
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        z += nm.exp(-1. * beta * e)
    return z

def factorized_q(beta, tar_seq, len_gen_seq, p_matrix, step_size=2., match_assured=False):  # Mean field value of boltzmann probability of match Q
    q = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        q_b = p_matrix[bases.index(tar_seq[i]), bases.index(tar_seq[i + 1])] * nm.exp(-1. * beta * pair_energy(pair, compl_seq(pair)))
        q = q * q_b
    q = q * (len_gen_seq - len(tar_seq) + 1) / 2  # Rescales Q so that it can be compared with the computed one

    if step_size == 1:
        q += factorized_q(beta, tar_seq, len_gen_seq, p_matrix, step_size=2, match_assured=False)
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        q += nm.exp(-1. * beta * e)
    return q

if __name__ == "__main__":
    A = 2
    if A == 1: # Plots selectivity comp vs selectivity math vs beta
        energy = rnd.uniform(-0.5, 0.5, (16, 16))  # binding energy matrix between all 16 base pairs
        nm.fill_diagonal(energy, 0)
        energy = energy + energy.T  # symmetrizing binding energy matrix
        p_bases = [0.25, 0.25, 0.25, 0.25]
        P = rnd.normal(1, 0.1, (4, 4))
        P = P / nm.sum(P)
        l = 4
        N_g = 10 ** 5

        x = nm.linspace(0, 10)
        t = random_seq(l, p_bases)
        g = random_seq_gen(N_g, P)

        Y = lambda x: q_perfect_match(x, t, g) / partition_function(x, t, g)
        y = lambda x: factorized_q(x, t, N_g, P) / factorized_z(x, t, N_g, P)
        plt.plot(x, Y(x), x, y(x))
        plt.show()

    if A == 2: # Calculates running average of Z comp vs Z math
        N_g = 10 ** 6
        beta = 1.623  # Physiological beta in kcal/mol
        energy = nm.loadtxt("energies1.txt")
        P = nm.loadtxt("prob_matrix.txt")
        g = random_seq_gen(N_g, P)
        t = "ATATATAT"
        z_math = factorized_z(beta, t, len(t), P, step_size=1)
        running_z = 0
        Y = []
        x = range(N_g - len(t))
        for i in x:
            z = partition_function(beta, t, g[i:i + len(t)], step_size=2)
            running_z += z
            y = running_z / (z_math * (i + 1))
            Y.append(y)
        #plt.ylim(0, 2)
        plt.title('t = (' + t + ')')
        plt.ylabel('Running average: Z comp / Z sym')
        plt.xlabel('Genome size')
        plt.plot(Y)
        plt.show()

    if A == 3: # basic comparison between the computated and the partitioned z's
        N_g = 10 ** 6
        beta = 1.623  # Physiological beta in kcal/mol
        energy = nm.loadtxt("energies1.txt")
        P = nm.loadtxt("prob_matrix.txt")
        g = random_seq_gen(N_g, P)
        t = "AAAAAAAA"
        step_size = 1
        z_math = factorized_z(beta, t, N_g, P, step_size=step_size)
        z_comp = partition_function(beta, t, g, step_size=step_size)
        print z_math/z_comp