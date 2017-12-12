from model3 import partition_function, q_perfect_match, factorized_z, factorized_q, random_seq
import numpy as nm
import model3
import matplotlib.pyplot as plt
beta = 1.

def log_error(num, denom): # handy function that returns the log of the fractional error
    return nm.abs(nm.log(num / denom))

def error(size_tar, size_gen):
    p_bases_tar = [0.25] * 4
    p_bases_gen = nm.loadtxt("prob_matrix.txt")
    tar = random_seq(size_tar, p_bases_tar)
    gen = random_seq(size_gen, p_bases_gen.sum(0))
    z_comp = partition_function(beta, tar_seq=tar, gen_seq=gen)
    q_comp = q_perfect_match(beta, tar_seq=tar, gen_seq=gen)
    z_math = factorized_z(beta, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen)
    q_math = factorized_q(beta, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen)
    return q_math * z_comp / (q_comp * z_math)

def selectivity_comp(size_tar, size_gen):
    p_bases_tar = [0.25] * 4
    p_bases_gen = nm.full((4, 4), 1. / 16.)
    tar = random_seq(size_tar, p_bases_tar)
    gen = random_seq(size_gen, p_bases_gen.sum(0))
    z_comp = partition_function(beta, tar_seq=tar, gen_seq=gen)
    q_comp = q_perfect_match(beta, tar_seq=tar, gen_seq=gen)
    # print q_comp
    return z_comp/q_comp

def selectivity_math(size_tar, size_gen):
    p_bases_tar = [0.25] * 4
    p_bases_gen = nm.full((4, 4), 1. / 16.)
    tar = random_seq(size_tar, p_bases_tar)
    z_math = factorized_z(beta, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen)
    q_math = factorized_q(beta, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen)
    return z_math/q_math

if __name__ == "__main__":
    A = 1
    if A == 1:
        model3.match_assured = True
        model3.energy = nm.loadtxt("energies1.txt")
        size_gen = 10**6
        sizes_tar = [4,6,8,10,12,14,16]
        trials = 30
        results = nm.zeros((len(sizes_tar), trials))
        for i, l in enumerate(sizes_tar):
            print l
            for j in range(trials):
                print j,
                results[i][j] = error(l, size_gen)
            print
        #print results
        nm.save("error_of_s_vs_l_10^6_bases_30_trials_4-16", results)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = nm.array(sizes_tar).repeat(trials)
        y = results.flatten()
        ax.scatter(x, y, s=10)
        ax.set_ylim(bottom=0, top=None)
        ax.xlabel("Targeting seq length (l)")
        ax.ylabel("S_math / S_comp")
        plt.show()
    if A == 2:
        size_gen = 10 ** 6
        sizes_tar = [4, 6, 8, 10, 12, 14, 16]
        trials = 10
        results = nm.zeros((len(sizes_tar), trials))
        results1= nm.zeros((len(sizes_tar), trials))
        for i, l in enumerate(sizes_tar):
            print l
            for j in range(trials):
                results[i][j] = selectivity_comp(l, size_gen)
        for i, l in enumerate(sizes_tar):
            print l
            for j in range(trials):
                results1[i][j] = selectivity_math(l, size_gen)
        # nm.save("selectivity_vs_l_10^5_bases_30_trials_16-20", results)
        fig = plt.figure()
        ax = fig.add_subplot(2, 1, 1)
        ax1 = fig.add_subplot(2, 1, 2)
        x = nm.array(sizes_tar).repeat(trials)
        y = results.flatten()
        y1 = results1.flatten()
        y2 = y1 / y # maths divided by comp
        ax.scatter(x, y, c='b', s=10)
        ax.scatter(x, y1, c='r', s=10)
        ax1.scatter(x, y2, c='g', s=10)
        #ax.set_ylim(bottom=0, top=None)
        plt.show()