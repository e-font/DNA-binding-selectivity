import numpy as nm
import numpy.random as rnd
from model3 import compl_seq

class Model3:
    bases = ['A', 'C', 'G', 'T']
    pairs = nm.array([[p + q for q in bases] for p in bases])

    def __init__(self, beta, match_assured=True):
        self.beta = beta
        self._energies = rnd.uniform(0, 0.5, (16,16))
        nm.fill_diagonal(self._energies, 0)
        self._energies = self._energies + self._energies.T
        self.match_assured = match_assured
        self.p_targeting = nm.full((4,), 1. / 4.)
        self.p_genome = nm.full((4, 4), 1. / 16.)
        self.t = None
        self.g = None
        self.matches = int(match_assured)
        self._z = 0

    def seq_set_check(self):
        if self.t is None or self.g is None:
            raise Exception("One or both of the sequences have not been initialized")

    def new_targeting(self, length):
        _t = self.Sequence()
        _t.seq = rnd.choice(Model3.bases, size=length, p=map(lambda x: float(x), self.p_targeting))
        self.t = _t
        self.matches = int(self.match_assured)
        self._z = 0

    def new_genome(self, length):
        _g = self.Sequence()
        _g.seq = rnd.choice(Model3.pairs, size=length, p=map(lambda x: float(x), self.p_genome))
        self.g = _g
        self.matches = int(self.match_assured)
        self._z = 0

    def self_energy(self, seq):
        e = 0.
        for p in seq.pairs():
            e += self._energies(p, compl_seq(p))
        return e

    class Sequence:
        def __init__(self):
            self._seq = []
            self.length = 0

        def pairs(self):
            i = 0
            while i < self.length:
                yield self._seq[i] + self.seq[i + 1]
                i += 2

        @property
        def seq(self):
            return self._seq

        @seq.setter
        def seq(self, sequence):
            self._seq = sequence
            self.length = len(sequence)

    def partition_function(self):  # actual partition function
        self.seq_set_check()
        _z = 0.
        for i in range(0, self.g.length - self.t.length, 2):
            e = 0.  # total energy for a given position of tar seq along gen seq
            matches = True
            for j in range(0, self.t.length, 2):  # stepping by two
                p1 = self.t.seq[j] + self.t.seq[j + 1]
                p2 = self.g.seq[i + j] + self.g.seq[i + j + 1]
                e += self._energies(p1, p2)
                if not p1 == p2:
                    matches = False
            if matches:
                self.matches += 1
            _z += nm.exp(-1. * self.beta * e)
        self._z = _z
        return _z

    def q_match(self):
        self.seq_set_check()
        if self.matches == int(self.match_assured) and self._z == 0:
            for i in range(0, self.g.length - self.t.length, 2):
                j = 0
                matches = True
                while matches and j < self.t.length / 2:  # stepping by two
                    p1 = self.t.seq[j] + self.t.seq[j + 1]
                    p2 = self.g.seq[i + j] + self.g.seq[i + j + 1]
                    if not p1 == p2:
                        matches = False
                    j += 1
                if matches:
                    self.matches += 1
        return self.matches * nm.exp(-1. * self.beta * self.self_energy(self.t))

m = Model3(beta=1)
m.new_targeting(5)
m.new_genome(1000)
print m.t.seq










