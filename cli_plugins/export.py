from cli_plugins.base.CLI_plugin import CLIPlugin

class ExportPlugin(CLIPlugin):

    def __init__(self, parser):
        super().__init__(parser, command_name="export", help_text = "Exporter le graphe d'une molécule à partir d'un ID Chebi")
        self.add_argument("--chebi_id", help_text="L'ID Chebi de la molécule", required=True)

    def execute(self, namespace, **kwargs):

        chebi_id = namespace.chebi_id

        print(f"Pas implémenté encore mais {chebi_id} a été reçu")