from model3 import partition_function, q_perfect_match, factorized_z, factorized_q, random_seq, random_seq_gen
import numpy as nm
import model3
import matplotlib.pyplot as plt
BETA = 1.623 #Physiological beta

def log_error(num, denom): # handy function that returns the log of the fractional error
    return nm.abs(nm.log(num / denom))

def error(tar, gen):
    steps = (len(gen) - len(tar) + 1)/2. #default step_size is 2
    z_comp = partition_function(BETA, tar_seq=tar, gen_seq=gen)
    #q_comp = q_perfect_match(BETA, tar_seq=tar, gen_seq=gen)
    z_math = factorized_z(BETA, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen) * steps
    #q_math = factorized_q(BETA, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen) * steps
    return z_comp / z_math

def selectivity_comp(tar, gen):
    tar = random_seq(size_tar, p_bases_tar)
    gen = random_seq(size_gen, p_bases_gen.sum(0))
    z_comp = partition_function(BETA, tar_seq=tar, gen_seq=gen)
    q_comp = q_perfect_match(BETA, tar_seq=tar, gen_seq=gen)
    # print q_comp
    return z_comp/q_comp

def selectivity_math(size_tar, size_gen):
    p_bases_tar = [0.25] * 4
    p_bases_gen = nm.full((4, 4), 1. / 16.)
    tar = random_seq(size_tar, p_bases_tar)
    z_math = factorized_z(BETA, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen)
    q_math = factorized_q(BETA, tar_seq=tar, len_gen_seq=size_gen, p_matrix=p_bases_gen)
    return z_math/q_math

if __name__ == "__main__":
    A = 2
    if A == 1:
        model3.match_assured = True
        model3.energy = nm.loadtxt("energies1.txt")
        size_gen = 10 ** 4
        sizes_tar = [4,6,8,10,12,14,16]
        p_bases_tar = [0.25] * 4
        p_bases_gen = nm.loadtxt("prob_matrix.txt")
        trials = 3
        results = nm.zeros((len(sizes_tar), trials))
        for i, l in enumerate(sizes_tar):
            print l
            for j in range(trials):
                #print j,
                tar = random_seq(l, p_bases_tar)
                gen = random_seq_gen(size_gen, p_bases_gen)
                results[i][j] = error(tar, gen)
            print
        #print results
        #nm.save("error_of_s_vs_l_10^6_bases_30_trials_4-16", results)
        #fig = plt.figure()
        #ax = fig.add_subplot(1, 1, 1)
        x = nm.array(sizes_tar).repeat(trials)
        y = results.flatten()
        plt.scatter(x, y, s=10)
        plt.ylim(bottom=0, top=None)
        plt.xlabel("Targeting seq length (l)")
        plt.ylabel("S_math / S_comp")
        plt.show()
    if A == 2:
        size_gen = 10 ** 5
        sizes_tar = [4, 6, 8, 10, 12, 14, 16]
        trials = 10
        results = nm.zeros((len(sizes_tar), trials))
        results1= nm.zeros((len(sizes_tar), trials))
        p_bases_tar = [0.25] * 4
        p_bases_gen = nm.loadtxt("prob_matrix.txt")
        model3.energy = nm.loadtxt("energies1.txt")
        for i, l in enumerate(sizes_tar):
            print l
            for j in range(trials):
                t = random_seq(l, p_bases_tar)
                g = random_seq_gen(size_gen, p_bases_gen)
                results[i][j] = partition_function(BETA, t, g)
                results1[i][j] = factorized_z(BETA, t, g, p_bases_gen) * ((size_gen - l + 1) / 2.)
        # nm.save("selectivity_vs_l_10^5_bases_30_trials_16-20", results)
        fig = plt.figure()
        ax = fig.add_subplot(2, 1, 1)
        ax1 = fig.add_subplot(2, 1, 2)
        x = nm.array(sizes_tar).repeat(trials)
        y = results.flatten()
        y1 = results1.flatten()
        y2 = y1 / y # maths divided by comp
        ax.scatter(x, nm.log(y), c='b', s=10)
        ax.scatter(x, nm.log(y1), c='r', s=10)
        ax1.scatter(x, y2, c='g', s=10)
        #ax.set_ylim(bottom=0, top=None)
        ax1.set_xlabel('Targeting seq length (l)')
        ax.set_ylabel('ln Z')
        ax1.set_ylabel('Z comp / Z math')
        fig.suptitle('Z comp vs Z vs l')
        plt.show()