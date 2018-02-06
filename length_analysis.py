'''
Script that investigates the accuracy of the mathematical model vs the computational model for different values of the
targeting sequence length l.
Imports from model3.py functions.
'''

from model3 import partition_function, q_perfect_match, factorized_z, factorized_q, random_seq, random_seq_gen
import numpy as nm
import model3
import matplotlib.pyplot as plt

# Variables
beta = 1.623  # Physiological beta
size_gen = 10 ** 5
sizes_tar = [4, 6, 8, 10, 12, 14, 16]  # Values of l to be investigated
trials = 10

# Creating useful arrays: results, probabilities, energy matrix
results_comp = nm.zeros((len(sizes_tar), trials))  # Will store the computed values of Z from the actual t and g
results_math = nm.zeros((len(sizes_tar), trials))  # Will store the modelled values of Z from the mathematical model
p_bases_tar = [0.25] * 4  # Uniform distribution for the bases in t
p_bases_gen = nm.loadtxt("prob_matrix.txt")
model3.energy = nm.loadtxt("energies1.txt")

# Calculating the values of z
for i, l in enumerate(sizes_tar):
    print l
    for j in range(trials):
        t = random_seq(l, p_bases_tar)  # Targeting sequence
        g = random_seq_gen(size_gen, p_bases_gen)  # Genome sequence
        results_comp[i][j] = partition_function(beta, t, g)
        results_math[i][j] = factorized_z(beta, t, g, p_bases_gen) * ((size_gen - l + 1) / 2.)  # Takes into account the step size

# Displays the results
fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)
ax1 = fig.add_subplot(2, 1, 2)
x = nm.array(sizes_tar).repeat(trials)
y = results_comp.flatten()
y1 = results_math.flatten()
y2 = y1 / y  # Z maths divided by Z comp - measure of error
ax.scatter(x, nm.log(y), c='b', s=10)
ax.scatter(x, nm.log(y1), c='r', s=10)
ax1.scatter(x, y2, c='g', s=10)
ax1.set_xlabel('Targeting seq length (l)')
ax.set_ylabel('ln Z')
ax1.set_ylabel('Z comp / Z math')
fig.suptitle('Z comp vs Z vs l')
plt.show()