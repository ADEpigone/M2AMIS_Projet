from cli_plugins.base.CLI_plugin import CLIPlugin

class IsomorphismTestPlugin(CLIPlugin):

    def __init__(self, parser):
        super().__init__(parser, command_name="iso", help_text = "Tester si deux molécules sont isomorphes")
        self.command_name = "iso"
        self.add_argument("--id1", help_text="L'ID Chebi de la molécule", required=True)
        self.add_argument("--id2", help_text="L'ID Chebi de la molécule", required=True)

    def execute(self, namespace, **kwargs):

        id1 = namespace.id1
        id2 = namespace.id2
        print(f"Pas implémenté encore mais les ids {id1} et {id2} ont été reçus")