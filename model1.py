"""
This is a program simulating the interaction between a CRISPR-like targeting DNA sequence, and a target DNA sequence.
The target is longer than the targeting one, and randomized with settable frequencies for the bases.
The program returns a value for the selectivity, which is defined as 1/(1/q - 1) where q is the energetic ratio between
a matching sequence and the sum of all sequences (the partition function Z).
"""
import numpy as nm
from numpy import random as rnd
import matplotlib.pyplot as plt


def hamiltonians(source_seq, target_seq, step=1):
    lst = []
    len_src = nm.size(source_seq)
    len_tar = nm.size(target_seq)
    for i in range(0, len_tar - len_src, step):
        m = 0
        for j in range(len_src):
            if target_seq[i + j] == source_seq[j]:
                m += 1
        lst.append((m, len_src - m))
    return lst


def partition_func(ham, e_match, e_no_match, b):
    z = sum(nm.exp(-1.*(i[0]*e_match+i[1]*e_no_match)*b) for i in ham)
    #print 'z: ', z,
    return z


def q_best_match(h, e_match, e_no_match, b):
    hams = nm.array(h)
    best = max(hams[:,0])
    num_of_bests = 0
    not_best = 0
    for item in hams:
        if item[0] == best:
            num_of_bests += 1
            not_best = item[1]

    q = num_of_bests * nm.exp(-1 * (e_match * best + e_no_match * not_best) * b)
    return q

def q_perfect_match(hams, l, e_match, b):
    num_of_match = 0
    for item in hams:
        if item[1] == 0:
            num_of_match += 1
    q = num_of_match * nm.exp(-1 * (e_match * l) * b)
    # print num_of_match
    return q


def selectivity(part_fun, q):
    return 1. / (part_fun / q - 1.)


def prob_tot(occurr_lst, bases_lst):
    return nm.prod(map(lambda i, j: float(i)**j, bases_lst[1], occurr_lst))


def high_beta_approx(b, e_match, e_no_match, occurr_lst, bases_lst):
    summ = sum(map(lambda i, j: i * (1. - float(j)) / float(j), occurr_lst, bases_lst[1]))
    return nm.exp((e_no_match - e_match)*b)/summ


def low_beta_approx(b, e_match, e_no_match, occurr_lst, bases_lst):
    p = prob_tot(occurr_lst, bases_lst)
    return p*(1-l*e_match*b)/(1-p-e_no_match*b*(sum(float(bases_lst[1,i])*occurr_lst[i] for i in range(4))-l*p))


def factorized_part_func(b, bases_lst, e_match, e_no_match, occurr_lst):
    return nm.prod(map(lambda p, o: (1. + float(p) * (nm.exp((e_no_match-e_match)*b) - 1.)) ** o, bases_lst[1], occurr_lst))


def factorized_q(b, bases_lst, e_match, e_no_match, occurr_lst):
    return nm.prod(map(lambda p, o: (float(p) * nm.exp((e_no_match-e_match)*b)) ** o, bases_lst[1], occurr_lst))


def factorized_selectivity(b, bases_lst, e_match, e_no_match, occurr_lst):
    z = factorized_part_func(b, bases_lst, e_match, e_no_match, occurr_lst)
    q = factorized_q(b, bases_lst, e_match, e_no_match, occurr_lst)
    return selectivity(z, q)

def correction_selectivity(s, p_match, N_g, l):
    return s * (1. + (p_match * (N_g - l)) ** -1.)

if __name__ == '__main__':
    BASES = nm.array([['A', 'T', 'C', 'G'], [0.2922, 0.2928, 0.2074, 0.2076]])
    #BASES = nm.array([['A', 'T', 'C', 'G'], [0.25, 0.25, 0.25, 0.25]])
    l = 10
    N_g = 100000
    E_MATCH = -1.
    E_NO_MATCH = 0.

    beta_min = 0
    beta_max = 1
    beta_num = 10
    x = nm.linspace(beta_min, beta_max, num=beta_num)
    source_sequence = rnd.choice(BASES[0], size=l, p=map(lambda x: float(x), BASES[1]))
    target_sequence = rnd.choice(BASES[0], size=N_g, p=map(lambda x: float(x), BASES[1]))
    target_sequence = nm.concatenate((source_sequence, target_sequence))

    occurrences = map(lambda b: sum(b == j for j in source_sequence), BASES[0])
    h = hamiltonians(source_sequence, target_sequence)

    #Y = nm.vectorize(lambda x: selectivity(partition_func(h, E_MATCH, E_NO_MATCH, x), q_best_match(h, E_MATCH, E_NO_MATCH, x)))(x)
    Y = nm.vectorize(
        lambda x: selectivity(partition_func(h, E_MATCH, E_NO_MATCH, x), q_perfect_match(h, l, E_MATCH, x)))(x)
    y = nm.vectorize(lambda x: factorized_selectivity(x, BASES, E_MATCH, E_NO_MATCH, occurrences))(x)

    # h = nm.vectorize(lambda x: high_beta_approx(x, E_MATCH, E_NO_MATCH, occurrences, BASES))(x)
    # l = nm.vectorize(lambda x: low_beta_approx(x, E_MATCH, E_NO_MATCH, occurrences, BASES))(x)
    # fig, ax = plt.subplots()
    #
    # ax.plot(x, Y/y, '-r', label='sym/factorized')
    # ax.plot(x, Y/l, '-g', label='sym/low-beta approx.')
    # ax.plot(x, Y/h, '-b', label='sym/high-beta approx.')
    # legend = ax.legend(loc='lower right', shadow=False)
    # plt.ylim(0., 1.2)
    # plt.title('Model 1: n=3, N=10^5, p=genome')
    # plt.show()
    p = prob_tot(occurrences, BASES)
    y1 = nm.vectorize(lambda x: correction_selectivity(x, p, N_g, l))(y)
    fig, ax = plt.subplots()
    ax.plot(x, Y, '-r', label='actual')
    ax.plot(x, y, '-g', label='without correction')
    ax.plot(x, y1, '-b', label='with correction')
    legend = ax.legend(loc='lower right', shadow=False)
    plt.show()
