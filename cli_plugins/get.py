from Chebi.CheBi2 import CheBi2
from cli_plugins.base.CLI_plugin import CLIPlugin
from utils import normalize_chebi_id

class GetPlugin(CLIPlugin):

    def __init__(self, parser, chebi_client: CheBi2 = None, **kwargs):
        super().__init__(
            parser,
            command_name="get",
            help_text="Récupérer une molécule depuis ChEBI et l'ajouter à chebi2",
        )
        self.add_argument("--chebi-id", help_text="ID ChEBI de la molécule", required=True)
        self.chebi_client = chebi_client

    def execute(self, namespace, **kwargs):
        if self.chebi_client is None:
            print("Client CheBi2 indisponible.")
            return

        chebi_id = normalize_chebi_id(namespace.chebi_id)

        print(f"Source sélectionnée: ChEBI ({chebi_id})")
        print("Récupération de la molécule...")
        mol_file = self.chebi_client.get_mol(chebi_id)
        if not mol_file:
            print(f"Aucune molécule trouvée pour {chebi_id}.")
            return

        print(f"Préparation insertion DB: chebi_id={chebi_id}")
        print("Mise à jour de la base chebi2...")
        self.chebi_client.update_table([(chebi_id, None, mol_file)])
        print(f"Molécule {chebi_id} récupérée depuis ChEBI et stockée dans chebi2.")