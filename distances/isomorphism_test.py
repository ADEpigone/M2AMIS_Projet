from distances.base.graph_distance import BaseGraphDistance
from graph import MoleculeGraph


class IsomorphismTest(BaseGraphDistance):
    def calculate_distance(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        """
        Test d'iso à mettre ici, retourne 0 si pas iso, 1 si iso ?
        On part du principe que iso = distance binaire, mais une distance !
        """
        if True:
            return 0.0
        else:
            return 1.0