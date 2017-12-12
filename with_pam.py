import numpy as nm
import numpy.random as rnd
import model1
import matplotlib.pyplot as plt

A = 2

E_MATCH = -1.
E_NO_MATCH = 0.
PAM = 'GG'

def partition_func(x, tar_seq, gen_seq):
    z = 0.
    len_pam = len(PAM)

    for i in range(len(gen_seq) - len(tar_seq)):
        if all(gen_seq[i + k] == PAM[k] for k in range(len_pam)):
            e = 0.
            for j in range(len(tar_seq)):
                if tar_seq[j] == gen_seq[i + len_pam + j]:
                    e += E_MATCH
                else:
                    e += E_NO_MATCH
            z += nm.exp(-1. * x * e)
    return z

def q_perfect_match(x, tar_seq, gen_seq):
    z = 0.

    len_pam = len(PAM)
    for i in range(len(gen_seq) - len(tar_seq)):
        if all(gen_seq[i + k] == PAM[k] for k in range(len_pam)):
            if all(gen_seq[i + len_pam + k] == tar_seq[k] for k in range(len(tar_seq))):
                z += nm.exp(-1. * x * E_MATCH * len(tar_seq))

    if z == 0:
        raise Warning('No matches!')
    return z

def pam_finder(gen_seq):
    len_pam = len(PAM)
    match_count = 0
    _space_counter = 0
    spaces_lst = []
    i = 0
    while i < (len(gen_seq) - len_pam + 1):
        if all(gen_seq[i + k] == PAM[k] for k in range(len_pam)):
            match_count += 1
            if _space_counter != 0:
                spaces_lst.append(i - len_pam - _space_counter)
            _space_counter = i
            i += 1
        i += 1

    return match_count, spaces_lst

def separ_distr(separation, p_pam):
    return 1 - (1 - p_pam) ** (separation)

if __name__ == '__main__':
    #BASES = nm.array([['A', 'T', 'C', 'G'], [0.2922, 0.2928, 0.2074, 0.2076]])
    BASES = nm.array([['A', 'T', 'C', 'G'], [0.25, 0.25, 0.25, 0.25]])
    l = 6
    N_g = 100000
    t = rnd.choice(BASES[0], size=l, p=map(lambda x: float(x), BASES[1]))
    g = rnd.choice(BASES[0], size=N_g, p=map(lambda x: float(x), BASES[1]))
    #g = 'AAAAGGAAAAATTAAACGGGAAAAGGAGGAAAACCAAGAAAGGAA'
    p_pam = 1./16.

    if A == 1:
        #print pam_finder(g)
        #print N_g * p_pam
        sep_list = nm.array(pam_finder(g)[1])
        #sep_list = sep_list/len(sep_list)
        y, x, patches = plt.hist(sep_list, bins=50, cumulative=True, color='b')
        #print x, y
        plt.show()
        y = nm.insert(y, 0, 0)

        plt.plot(x, y/y[-1], x, (1 - (1 - 1./16.) ** x))
        plt.show()

    if A == 2:
        h = model1.hamiltonians(t, g)
        Y = nm.vectorize(lambda x: model1.selectivity(partition_func(x, t, g), q_perfect_match(x, t, g)))
        y = nm.vectorize(lambda x: model1.selectivity(model1.partition_func(h, E_MATCH, E_NO_MATCH, x), model1.q_perfect_match(h, l, E_MATCH, x)))
        x = nm.linspace(0, 10)
        #print Y(10), partition_func(10, t, g), q_perfect_match(10, t, g)
        plt.plot(x, Y(x)/y(x))
        plt.show()