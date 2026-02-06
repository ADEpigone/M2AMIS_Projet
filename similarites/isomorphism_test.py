from similarites.base.graph_similarity import BaseGraphSimilarity
from graph import MoleculeGraph


class IsomorphismTest(BaseGraphSimilarity):
    def calculate_similarity(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        """
        Test d'iso à mettre ici, retourne 0 si pas iso, 1 si iso ?
        On part du principe que iso = distance binaire, mais une distance !
        """
        if True:
            return 0.0
        else:
            return 1.0