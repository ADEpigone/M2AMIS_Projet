from Chebi.CheBi2 import CheBi2
from cli_plugins.base.CLI_plugin import CLIPlugin

from Iso.iso_test import IsoTest
from graph import MoleculeGraph
from utils import check_none

class IsomorphismTestPlugin(CLIPlugin):

    def __init__(self, parser, chebi_client : CheBi2, **kwargs):
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

        if check_none(mol1, id1) or check_none(mol2, id2):
            return
        
        chebi_id1 = id1 if str(id1).upper().startswith("CHEBI:") else f"CHEBI:{id1}"
        chebi_id2 = id2 if str(id2).upper().startswith("CHEBI:") else f"CHEBI:{id2}"
        print(f"Molécule 1: https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id1}")
        print(f"Molécule 2: https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id2}")

        g1 = MoleculeGraph.from_moltext(mol1)
        g2 = MoleculeGraph.from_moltext(mol2)

        if self.test.calculate_similarity(g1, g2):
            print(f"Les molécules {id1} et {id2} sont isomorphes.")
        else:
            print(f"Les molécules {id1} et {id2} ne sont pas isomorphes.")