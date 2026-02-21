from Chebi.CheBi import Chebi
from cli_plugins.base.CLI_plugin import CLIPlugin
from graph import MoleculeGraph

from similarites.cwl_kernel import CWLKernel
from similarites.builtin_similarity import BuiltinSimilarity
from utils import print_ids

class ComparisonPlugin(CLIPlugin):

    def __init__(self, parser, chebi_client : Chebi, **kwargs):
        super().__init__(parser, command_name="comp", help_text = "Comparer la similarité de deux molécules")
        self.add_argument("--id1", help_text="L'ID Chebi de la molécule", required=True)
        self.add_argument("--id2", help_text="L'ID Chebi de la molécule", required=True)
        self.add_argument("--fingerprint", help_text="Le type de fingerprint à utiliser pour la comparaison", required=True)
        self.add_argument("--method", help_text="La méthode de comparaison à utiliser (tanimoto, etc.)", required=True)
        self.chebi_client = chebi_client

    def execute(self, namespace, **kwargs):

        id1 = namespace.id1
        id2 = namespace.id2
        fingerprint = namespace.fingerprint.lower()
        method = namespace.method

        mol1 = self.chebi_client.get_mol(id1)
        mol2 = self.chebi_client.get_mol(id2)

        g1 = MoleculeGraph.from_moltext(mol1)
        g2 = MoleculeGraph.from_moltext(mol2)

        if fingerprint == "cwl":
            sim_m = CWLKernel(similarity=method)
        else:
            sim_m = BuiltinSimilarity(fingerprint=fingerprint, similarity=method)
        print_ids(id1, id2)
        print("Méthode de fingerprint : " + fingerprint.capitalize())
        print("Méthode de similarité : " + method.capitalize())
        print(f"La similarité entre les molécules {id1} et {id2} est : {sim_m.calculate_similarity(g1, g2)}")