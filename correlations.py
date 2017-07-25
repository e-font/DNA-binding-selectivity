import numpy as nm
import model1
import matplotlib.pyplot as plt
import re

A = 1
bases = ['A', 'T', 'C', 'G']
pairs = nm.array([[p + q for q in bases] for p in bases])

def factorized_z(x, tar_seq, p_matrix):
    z = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        z_bp = 0. #z for the base pair
        for j, row in enumerate(pairs):
            for k, p in enumerate(row):
                z_bp += p_matrix[j, k] * nm.exp(x * sum(p[i] == pair[i] for i in (0, 1)))
        z = z * z_bp
    return z

def factorized_q(x, tar_seq, p_matrix):
    z = 1.
    for i in range(0, len(tar_seq), 2):
        pair = tar_seq[i] + tar_seq[i + 1]
        for j, row in enumerate(pairs):
            for k, p in enumerate(row):
                if pair == p:
                    z = z*p_matrix[j, k]
    return z * nm.exp(x * len(tar_seq))

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

if __name__ == '__main__':
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
        x = nm.linspace(0, 2)

        #print P, seq_analyse(g)
        z_a = lambda i: z(i, d, l, 0, 0)
        z_b = lambda i: z(i, d, 0, l, 0)
        z_c = lambda i: z(i, d, 0, 0, l)

        q_a = lambda i: q(i, d, l, 0, 0)
        q_b = lambda i: q(i, d, 0, l, 0)
        q_c = lambda i: q(i, d, 0, 0, l)

        h_a1 = model1.hamiltonians(list(a), list(g))
        h_a2 = model1.hamiltonians(list(a), list(g), step=2)
        h_b = model1.hamiltonians(list(b), list(g))
        h_c = model1.hamiltonians(list(c), list(g))

        # print model1.partition_func(h_a1, e_match, e_no_match, 10), model1.partition_func(h_a2, e_match, e_no_match, 10),
        # print model1.q_perfect_match(h_a1, l*2, e_match, 10), model1.q_perfect_match(h_a2, l*2, e_match, 10)
        # print model1.partition_func(h_a1, e_match, e_no_match, 10) / model1.q_perfect_match(h_a1, l*2, e_match, 10), model1.partition_func(h_b, e_match, e_no_match, 10) / model1.q_perfect_match(h_b, l*2, e_match, 10)
        # print (z_a(10) + z_c(10)) / (q_a(10) + q_c(10)), (z_b(i) + z_c(i))/(q_b(i) + q_c(i))

        Y_a = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_a1, e_match, e_no_match, i),
                                                        model1.q_perfect_match(h_a1, l*2, e_match, i)))
        Y_b = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_b, e_match, e_no_match, i),
                                                        model1.q_perfect_match(h_b, l*2, e_match, i)))
        Y_c = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_c, e_match, e_no_match, i),
                                                        model1.q_perfect_match(h_c, l*2, e_match, i)))

        y_a = nm.vectorize(lambda i: model1.selectivity(z_a(i) + z_c(i), q_a(i) + q_c(i)))
        y_b = nm.vectorize(lambda i: model1.selectivity(z_b(i) + z_c(i), q_b(i) + q_c(i)))
        y_c = nm.vectorize(lambda i: model1.selectivity(z_c(i) * 2, q_c(i) * 2))

        fact_y_a = nm.vectorize(lambda i: model1.selectivity(factorized_z(i, a, P) + factorized_z(i, c, P), factorized_q(i, a, P) + factorized_q(i, c, P)))
        fact_y_b = nm.vectorize(lambda i: model1.selectivity(factorized_z(i, b, P) + factorized_z(i, c, P), factorized_q(i, b, P) + factorized_q(i, c, P)))
        fact_y_c = nm.vectorize(lambda i: model1.selectivity(factorized_z(i, c, P) * 2, factorized_q(i, c, P) * 2))

        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        f.suptitle('n_g = 10^6, d = -0.5')
        ax1.set_title('All matching pairs')
        ax1.plot(x, Y_a(x), x, y_a(x), x, fact_y_a(x))
        ax2.set_title('One base per pair')
        ax2.plot(x, Y_b(x), x, y_b(x), x, fact_y_b(x))
        ax3.set_title('No matching bases')
        ax3.plot(x, Y_c(x), x, y_c(x), x, fact_y_c(x))
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

    if A == 6:
        anomaly = (2, 3)
        e_match = -1.
        e_no_match = 0.
        l = 3  # Number of pairs in t

        a = pairs[anomaly] * l
        b = 'CA' * l
        c = 'AA' * l

        filename = 'hs_ref_GRCh38.p7_chr5.txt'
        with open(filename, 'r') as file:
            s = file.read()
            g = s[0:1000000]
        d = -0.86
        x = nm.linspace(0, 10)

        #print P, seq_analyse(g)
        z_a = lambda i: z(i, d, l, 0, 0)
        z_b = lambda i: z(i, d, 0, l, 0)
        z_c = lambda i: z(i, d, 0, 0, l)

        q_a = lambda i: q(i, d, l, 0, 0)
        q_b = lambda i: q(i, d, 0, l, 0)
        q_c = lambda i: q(i, d, 0, 0, l)

        h_a1 = model1.hamiltonians(list(a), list(g))
        h_a2 = model1.hamiltonians(list(a), list(g), step=2)
        h_b = model1.hamiltonians(list(b), list(g))
        h_c = model1.hamiltonians(list(c), list(g))

        # print model1.partition_func(h_a1, e_match, e_no_match, 10), model1.partition_func(h_a2, e_match, e_no_match, 10),
        # print model1.q_perfect_match(h_a1, l*2, e_match, 10), model1.q_perfect_match(h_a2, l*2, e_match, 10)
        # print model1.partition_func(h_a1, e_match, e_no_match, 10) / model1.q_perfect_match(h_a1, l*2, e_match, 10), model1.partition_func(h_b, e_match, e_no_match, 10) / model1.q_perfect_match(h_b, l*2, e_match, 10)
        # print (z_a(10) + z_c(10)) / (q_a(10) + q_c(10)), (z_b(i) + z_c(i))/(q_b(i) + q_c(i))

        Y_a = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_a1, e_match, e_no_match, i),
                                                        model1.q_perfect_match(h_a1, l*2, e_match, i)))
        Y_b = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_b, e_match, e_no_match, i),
                                                        model1.q_perfect_match(h_b, l*2, e_match, i)))
        Y_c = nm.vectorize(lambda i: model1.selectivity(model1.partition_func(h_c, e_match, e_no_match, i),
                                                        model1.q_perfect_match(h_c, l*2, e_match, i)))

        y_a = nm.vectorize(lambda i: model1.selectivity(z_a(i) + z_c(i), q_a(i) + q_c(i)))
        y_b = nm.vectorize(lambda i: model1.selectivity(z_b(i) + z_c(i), q_b(i) + q_c(i)))
        y_c = nm.vectorize(lambda i: model1.selectivity(z_c(i) * 2, q_c(i) * 2))

        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        f.suptitle('Hum genome seq: _g = 10^6, d = -0.86')
        ax1.set_title('All matching pairs')
        ax1.plot(x, Y_a(x), x, y_a(x))
        ax2.set_title('One base per pair')
        ax2.plot(x, Y_b(x), x, y_b(x))
        ax3.set_title('No matching bases')
        ax3.plot(x, Y_c(x), x, y_c(x))
        plt.show()

    if A == 7: # Just to check if partitioned z works
        anomaly = (2, 3) # CG
        l = 3  # Number of pairs in t
        a = pairs[anomaly] * l
        b = 'CA' * l
        c = 'AA' * l

        P = nm.full((4, 4), 1 / 16.)
        d = -0.5
        P += (-1. * d) / (15. * 16.)
        P[anomaly] += d / (15. * 16.) + d / 16.

        x = nm.linspace(0, 2)

        fact_y_a = nm.vectorize(lambda i: model1.selectivity(factorized_z(i, a, P), factorized_q(i, a, P)))
        fact_y_b = nm.vectorize(lambda i: model1.selectivity(factorized_z(i, b, P), factorized_q(i, b, P)))
        fact_y_c = nm.vectorize(lambda i: model1.selectivity(factorized_z(i, c, P), factorized_q(i, c, P)))

        z_a = lambda i: z(i, d, l, 0, 0)
        z_b = lambda i: z(i, d, 0, l, 0)
        z_c = lambda i: z(i, d, 0, 0, l)

        q_a = lambda i: q(i, d, l, 0, 0)
        q_b = lambda i: q(i, d, 0, l, 0)
        q_c = lambda i: q(i, d, 0, 0, l)

        y_a = nm.vectorize(lambda i: model1.selectivity(z_a(i), q_a(i)))
        y_b = nm.vectorize(lambda i: model1.selectivity(z_b(i), q_b(i)))
        y_c = nm.vectorize(lambda i: model1.selectivity(z_c(i), q_c(i)))

        # f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        # f.suptitle('n_g = 10^6, d = -0.5')
        # ax1.set_title('All matching pairs')
        # ax1.plot(x, fact_y_a(x), '-r', x, y_a(x), '-b')
        # ax2.set_title('One base per pair')
        # ax2.plot(x, fact_y_b(x), '-r', x, y_b(x), '-b')
        # ax3.set_title('No matching bases')
        # ax3.plot(x, fact_y_c(x), '-r', x, y_c(x), '-b')
        # plt.show()
        print factorized_z(2, a, P)
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        f.suptitle('n_g = 10^6, d = -0.5')
        ax1.set_title('All matching pairs')
        ax1.plot(x, factorized_z(x, a, P), '-r', x, z_a(x), '-b')
        ax2.set_title('One base per pair')
        ax2.plot(x, factorized_z(x, b, P), '-r', x, z_b(x), '-b')
        ax3.set_title('No matching bases')
        ax3.plot(x, factorized_z(x, c, P), '-r', x, z_c(x), '-b')
        plt.show()