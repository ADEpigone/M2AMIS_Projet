import rdkit
from similarites.base.graph_similarity import BaseGraphSimilarity
from abc import ABC, abstractmethod
from graph import MoleculeGraph
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

class BaseSimilarityFromFingerprint(BaseGraphSimilarity, ABC):

    def __init__(self, similarity_function=DataStructs.TanimotoSimilarity):
        self.similarity_function = similarity_function
    
    @abstractmethod
    def calculate_fingerprint(self, mol: MoleculeGraph) -> rdkit.DataStructs:
        pass

    def calculate_similarity(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        return self.similarity_function(self.calculate_fingerprint(g1), self.calculate_fingerprint(g2))
