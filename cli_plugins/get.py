
#Ici on est load dans le dossier d'avant donc il faut se replacer dans plugins
from cli_plugins.base.CLI_plugin import CLIPlugin
#Mais l'avantage c'est qu'on peut load utils super facilement donc GOATED
from utils import get_mol_file

class GetPlugin(CLIPlugin):

    def __init__(self, parser):
        super().__init__(parser, command_name="get", help_text = "Récupère le fichier mol d'une molécule à partir de son ID Chebi")
        self.add_argument("--chebi_id", help_text="L'ID Chebi de la molécule", required=True)

    def execute(self, namespace, **kwargs):

        chebi_id = namespace.chebi_id

        molfile = get_mol_file(chebi_id)
        if molfile:
            print(f"Fichier mol pour l'ID Chebi {chebi_id} :")
            print(molfile)
            #on return pas jsp quoi faire encore
        else:
            print(f"Aucun fichier mol trouvé pour l'ID Chebi {chebi_id}.")