import numpy as nm
import model1
import matplotlib.pyplot as plt
import re

A = 1
bases = ['A', 'T', 'C', 'G']
pairs = nm.array([[p + q for q in bases] for p in bases])

def z_0(x, pairs):
    return (0.25 * (nm.exp(x) + 3)) ** (pairs * 2.)


def q_0(x, pairs):
    return (0.25 * (nm.exp(x))) ** (pairs * 2)


def z(x, d, n_2, n_1, n_0):
    return z_0(x, n_2 + n_1 + n_0) * (1. + (d/15.) * (16. * nm.exp(2*x) / ((nm.exp(x) + 3) ** 2) - 1)) ** n_2 \
           * (1. + (d/15.) * (16. * nm.exp(x) / ((nm.exp(x) + 3) ** 2) - 1)) ** n_1 \
            * (1. + (d/15.) * (16. / ((nm.exp(x) + 3) ** 2) - 1)) ** n_0


def q(x, d, n_2, n_1, n_0):
    return q_0(x, n_2 + n_1 + n_0) * (1. + d) ** n_2 * (1. - (d/15.)) ** (n_1 + n_0)


def random_pair(P):
    i = nm.random.choice(range(4), 1, p=nm.sum(P, axis=1))[0]
    j = nm.random.choice(range(4), 1, p=P[i] / nm.sum(P[i]))[0]
    return pairs[i, j]


def seq_analyse(seq):
    prob_matrix = nm.zeros(pairs.shape)
    for i in range(0, len(seq), 2):
        pair = seq[i] + seq[i + 1]
        for j, row in enumerate(pairs):
            for k, p in enumerate(row):
                if pair == p:
                    prob_matrix[j, k] += 1
                    break
    prob_matrix = prob_matrix / float(prob_matrix.sum())
    return prob_matrix

if A == 1:
    anomaly = (2, 3)
    e_match = -1.
    e_no_match = 0.
    l = 3  # Number of pairs in t

    a = pairs[anomaly] * l
    b = 'CA' * l
    c = 'AA' * l

    N = 1000000 #Number of bases in g
    P = nm.full((4, 4), 1/16.)
    d = -0.5
    P += (-1. * d) / (15. * 16.)
    P[anomaly] += d / (15. * 16.) + d / 16. #CG
    g = ''
    for i in range(N / 2):
        g += random_pair(P)

    #print P, seq_analyse(g)
    z_a = lambda i: z(i, d, l, 0, 0)
    z_b = lambda i: z(i, d, 0, l, 0)
    z_c = lambda i: z(i, d, 0, 0, l)

    h_a = model1.hamiltonians(list(a), list(g), step=2)
    h_b = model1.hamiltonians(list(b), list(g), step=2)
    h_c = model1.hamiltonians(list(c), list(g), step=2)

    Y_a = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_a, e_match, e_no_match, i),
                                                    model1.q_perfect_match(h_a, l*2, e_match, i)))
    Y_b = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_b, e_match, e_no_match, i),
                                                    model1.q_perfect_match(h_b, l*2, e_match, i)))
    Y_c = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_c, e_match, e_no_match, i),
                                                    model1.q_perfect_match(h_c, l*2, e_match, i)))

    y_a = nm.vectorize(lambda i: model1.selectivity(z_a(i), q(i, d, l, 0, 0)))
    y_b = nm.vectorize(lambda i: model1.selectivity(z_b(i), q(i, d, 0, l, 0)))
    y_c = nm.vectorize(lambda i: model1.selectivity(z_c(i), q(i, d, 0, 0, l)))

    x = nm.linspace(0, 10)

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    f.suptitle('n_g = 10^6, d = -0.5')
    ax1.set_title('All matching pairs')
    ax1.plot(x, Y_a(x), x, y_a(x))
    ax2.set_title('One base per pair')
    ax2.plot(x, Y_b(x), x, y_b(x))
    ax3.set_title('No matching bases')
    ax3.plot(x, Y_c(x), x, y_c(x))
    plt.show()

if A == 3:
    path = 'C:\Users\Enrico\Downloads\\'
    infilename = 'hs_ref_GRCh38.p7_chr5.mfa'
    outfilename = 'hs_ref_GRCh38.p7_chr5.txt'
    with open(infilename, 'rb') as infile, open(outfilename, 'w') as outfile:
        lines = [line.upper().replace('N', '').replace('\n', '') if not re.search('ref', line) else '' for line in infile]
        outfile.writelines(lines)

if A == 4:
    filename = 'hs_ref_GRCh38.p7_chr5.txt'
    outfilename = 'prob_matrix.txt'
    with open(filename, 'r') as file, open(outfilename, 'w') as outfile:
        s = file.read()
        m = seq_analyse(s)
        outfile.write(str(m))

if A == 5:
    filename = 'prob_matrix.txt'
    m = nm.loadtxt(filename)
    print m
