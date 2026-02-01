from abc import ABC, abstractmethod
from typing import Any


class Distance(ABC):
    """
    Interface pour toutes les distances entre molécules.
    """

    @abstractmethod
    def calculer_distance(self, mol1: Any, mol2: Any) -> float:
        """
        Calcule une distance (ou dissimilarité) entre deux molécules.
        Plus la valeur est faible, plus les molécules sont similaires.
        """
        pass
