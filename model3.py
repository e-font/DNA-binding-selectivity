import numpy as nm
import numpy.random as rnd
import model1
import correlations
import matplotlib.pyplot as plt

bases = ['A', 'C', 'G', 'T']
pairs = nm.array([[p + q for q in bases] for p in bases])  # matrix of base pairs
energy = rnd.normal(0.5, 0.05, (16, 16))  # binding energy matrix between all 16 base pairs
energy = rnd.uniform(-0.5, 0.5, (16, 16))
nm.fill_diagonal(energy, 0)
energy = energy + energy.T  # symmetrizing binding energy matrix
# print energy
match_assured = False

def compl_seq(seq):
    c = ''
    for b in seq:
        i = int(-1. * bases.index(b)) - 1
        c += bases[i]
    return c

def pair_energy(pair1, pair2):  # returns binding energy given two strings (= base pairs)
    def pair_index(pair):
        return bases.index(pair[0]) * 4 + bases.index(pair[1])
    return energy[pair_index(pair1), pair_index(pair2)]

def partition_function(x, tar_seq, gen_seq):  # actual partition function
    z = 0.
    i = 0
    while i <= len(gen_seq) - len(tar_seq):
        #print 'IN'
        e = 0. # total energy for a given position of tar seq along gen seq
        for j in range(0, len(tar_seq), 2): # stepping by base pairs
            p1 = tar_seq[j] + tar_seq[j + 1]
            p2 = gen_seq[i + j] + gen_seq[i + j + 1]
            e += pair_energy(p1, p2)
        z += nm.exp(-1. * x * e)
        i += 2
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        z += nm.exp(-1. * x * e)
    return z

def q_perfect_match(x, tar_seq, gen_seq):  # actual number of perfect matches, which as e=0 is q
    q = 0.
    for i in range(0, len(gen_seq) - len(tar_seq), 2):
        if all(tar_seq[j] == gen_seq[j + i] for j in range(len(tar_seq))):
            e = 0.
            for j in range(0, len(tar_seq), 2):  # stepping by two
                p1 = tar_seq[j] + tar_seq[j + 1]
                p2 = gen_seq[i + j] + gen_seq[i + j + 1]
                e += pair_energy(p1, compl_seq(p2))
            q += nm.exp(-1 * x * e)
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        q += nm.exp(-1. * x * e)
    return q

def factorized_z(x, tar_seq, len_gen_seq, p_matrix):  # model relying on factorized equation for Z
    z = 1.
    for i in range(0, len(tar_seq), 2): #running over each pair in t
        pair = tar_seq[i] + tar_seq[i + 1]
        z_b = 0. #Z for the base pair
        for j, row in enumerate(pairs): # running over all 16 possible base pairs
            for k, p in enumerate(row):
                z_b += p_matrix[j, k] * nm.exp(-1. * x * pair_energy(pair, p))
        z = z * z_b
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        z += nm.exp(-1. * x * e) / (len_gen_seq - len(tar_seq) + 1)
    return z

def factorized_q(x, tar_seq, len_gen_seq, p_matrix):  # model relying on factorized equation for Q_match
    q = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        z_b = 0.
        for j, row in enumerate(pairs):
            for k, p in enumerate(row):
                if pair == p:
                    z_b = p_matrix[j, k] * nm.exp(-1. * x * pair_energy(p, compl_seq(p)))
                    break
        q = q * z_b
    if match_assured:
        e = 0.
        for i in range(0, len(tar_seq), 2):
            p = tar_seq[i] + tar_seq[i + 1]
            e += pair_energy(p, compl_seq(p))
        q += nm.exp(-1. * x * e) / (len_gen_seq - len(tar_seq) + 1)
    return q

def random_seq(size, p_bases):
    return rnd.choice(bases, size, p=map(lambda x: float(x), p_bases))

if __name__ == "__main__":
    A = 3
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
        y = lambda x: model1.selectivity(factorized_z(x, t[1:0:1], N_g, P) + factorized_z(x, t, N_g, P), factorized_q(x, t[1:0:1], N_g, P) + factorized_q(x, t, N_g, P))
        #factorized_q(t, P)
        #factorized_z(2, t, P)
        #print partition_function(10,t,g), q_perfect_match(10, t, g)
        #print factorized_z(x, t[1:0:1], P), factorized_z(x, t, P), factorized_q(t[1:0:1], P), factorized_q(t, P)
        plt.plot(x, Y(x))
        plt.show()
    if A == 2:
        print compl_pair('AT')

    if A == 3:
        N_g = 10 ** 5
        beta = 1.623 #beta in kcal/mol at 37C
        energy = nm.rot90(nm.loadtxt("energies1.txt"))
        g = ''
        #P = nm.loadtxt("prob_matrix.txt")
        #P = rnd.normal(1, 0.1, (4, 4))
        #P = P / nm.sum(P)
        P = nm.full((4,4), 1./16.)
        for i in range(N_g / 2):
            g += correlations.random_pair(P)
        t = "CGCGCGCG"
        z_math = factorized_z(beta, t, len(t), P)
        #print z_math
        running_z = 0
        Y = []
        x = range(N_g - len(t))
        for i in x:
            z = partition_function(beta, t, g[i:i + len(t)])
            #print z
            running_z += z
            #print g[i:i + len(t)]
            y = running_z / (z_math * (i + 1))
            Y.append(y)
            #print y
        print Y[-1]
        print running_z , (z_math * (N_g - len(t)))
        print correlations.seq_analyse(g) - P
        plt.ylim(0, 2)
        plt.ylabel('Running average: Z comp / Z sym')
        plt.xlabel('Genome size')
        plt.plot(Y)
        plt.show()


