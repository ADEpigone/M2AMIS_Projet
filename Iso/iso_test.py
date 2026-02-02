from Iso.vec_to_vc import *
from Iso.nauty_wrapper import *

class iso_test:

    def __init__(self, g1, g2):
        self.graph1 = g1
        self.graph2 = g2

    def are_isomorphic(self, directed):

        pg1 = to_pynauty(to_vc(self.graph1), directed)
        pg2 = to_pynauty(to_vc(self.graph2), directed)

        return "Are isomorphs ? : ", pynauty.certificate(pg1) == pynauty.certificate(pg2)


if __name__ == "__main__":
    g1 = MoleculeGraph.from_molfile("CHEBI_187524.mol")
    g2 = MoleculeGraph.from_molfile("CHEBI_187524.mol")

    test = iso_test(g1, g2)
    print(test.are_isomorphic(directed=False))