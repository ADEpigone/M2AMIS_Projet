from Iso.vec_to_vc import *
from Iso.nauty_wrapper import *
from similarites.base.graph_similarity import BaseGraphSimilarity

class IsoTest(BaseGraphSimilarity):
    def calculate_similarity(self, g1 : MoleculeGraph, g2: MoleculeGraph) -> float:
        if self.are_isomorphic(g1, g2, directed=False):
            return 1.0
        return 0.0

    def are_isomorphic(self,g1 : MoleculeGraph, g2: MoleculeGraph, directed = False) -> bool:

        pg1 = to_pynauty(to_vc(g1), directed)
        pg2 = to_pynauty(to_vc(g2), directed)

        return pynauty.certificate(pg1) == pynauty.certificate(pg2)


if __name__ == "__main__":
    g1 = MoleculeGraph.from_molfile("CHEBI_187524.mol")
    g2 = MoleculeGraph.from_molfile("CHEBI_187524.mol")

    test = IsoTest()
    print(test.calculate_similarity(g1, g2))