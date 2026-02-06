from abc import ABC, abstractmethod
from typing import Any
from graph import *

class BaseGraphSimilarity(ABC):
    @abstractmethod
    def calculate_similarity(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        pass
