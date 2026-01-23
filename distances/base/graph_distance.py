from abc import ABC, abstractmethod
from typing import Any
from graph import *

class BaseGraphDistance(ABC):
    @abstractmethod
    def calculate_distance(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        pass
