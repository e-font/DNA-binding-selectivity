import numpy as nm
import model1
import matplotlib.pyplot as plt
import re

A = 4
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


def random_pair():
    i = nm.random.choice(range(4), 1, p=nm.sum(P, axis=1))[0]
    j = nm.random.choice(range(4), 1, p=P[i] / nm.sum(P[i]))[0]
    return pairs[i, j]

def seq_analyse(seq):
    prob_matrix = nm.zeros(pairs.shape)
    for i in range(0, len(seq), 2):
        print i
        pair = seq[i] + seq[i + 1]
        prob_matrix[nm.where(pair == pairs)] += 2./len(seq)
    print prob_matrix

if A == 1:
    n = 5.
    n_2 = 5.
    n_1 = 0.
    n_0 = n - n_2 - n_1

    d = -0.8
    w_0 = lambda i: z_0(i, n)
    w = lambda i: z(x, d, n_2, n_1, n_0)
    x = nm.linspace(0, 0.2)
    y_0 = nm.vectorize(lambda i: model1.selectivity(w_0(i), q_0(i, n)))
    y = nm.vectorize(lambda i: model1.selectivity(w(i), q(i, d, n_2, n_1, n_0)))

    plt.plot(x, w_0(x), x, w(x))
    plt.show()

if A == 2:
    P = nm.full((4, 4), 1/16.)
    P[3, 2] -= 0.01
    P[2, 3] += 0.01
    s = ''
    for i in range(100000):
        s += random_pair()
    print P
    seq_analyse(s)

if A == 3:
    path = 'C:\Users\Enrico\Downloads\\'
    infilename = 'hs_ref_GRCh38.p7_chr5.mfa'
    outfilename = 'hs_ref_GRCh38.p7_chr5.txt'
    with open(infilename, 'rb') as infile, open(outfilename, 'w') as outfile:
        lines = [line.upper().replace('N', '').replace('\n', '') if not re.search('ref', line) else '' for line in infile]
        outfile.writelines(lines)

if A == 4:
    s = ''
    filename = 'hs_ref_GRCh38.p7_chr5.txt'
    with open(filename, 'r') as file:
        s = file.read()
    seq_analyse(s)
