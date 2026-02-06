from Chebi.CheBi import Chebi
from cli_plugins.base.CLI_plugin import CLIPlugin

from Iso.iso_test import IsoTest
from graph import MoleculeGraph


class IsomorphismTestPlugin(CLIPlugin):

    def __init__(self, parser, chebi_client : Chebi, **kwargs):
        super().__init__(parser, command_name="iso", help_text = "Tester si deux molécules sont isomorphes")
        self.add_argument("--id1", help_text="L'ID Chebi de la molécule", required=True)
        self.add_argument("--id2", help_text="L'ID Chebi de la molécule", required=True)

        self.chebi_client = chebi_client
        self.test = IsoTest()

    def execute(self, namespace, **kwargs):

        id1 = namespace.id1
        id2 = namespace.id2

        mol1 = self.chebi_client.get_mol(id1)
        mol2 = self.chebi_client.get_mol(id2)

        g1 = MoleculeGraph.from_moltext(mol1)
        g2 = MoleculeGraph.from_moltext(mol2)

        print(self.test.calculate_similarity(g1, g2))