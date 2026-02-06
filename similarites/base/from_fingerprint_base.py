from collections import Counter
import rdkit
from similarites.base.graph_similarity import BaseGraphSimilarity
from abc import ABC, abstractmethod
from graph import MoleculeGraph
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

class BaseSimilarityFromFingerprint(BaseGraphSimilarity, ABC):

    def __init__(self, similarity_function=DataStructs.TanimotoSimilarity, n_bits=2048):
        self.similarity_function = similarity_function
        self.n_bits = n_bits
    
    @abstractmethod
    def calculate_fingerprint(self, mol: MoleculeGraph) -> rdkit.DataStructs:
        pass

    def _counter_to_fingerprint(self, counter: Counter) -> DataStructs:
        """
        Convertit un Counter avec des clés string en une ExplicitBitVect.
        Chaque clé est hachée pour obtenir un index de bit dans [0, n_bits).
        """
        bit_vect = DataStructs.ExplicitBitVect(self.n_bits)
        for key, count in counter.items():
            if count > 0:
                index = hash(key) % self.n_bits
                bit_vect.SetBit(index)
        return bit_vect

    def calculate_similarity(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        return self.similarity_function(self.calculate_fingerprint(g1), self.calculate_fingerprint(g2))
