
class CLIPlugin:

    def __init__(self, parser, command_name, help_text="Pas d'aide disponible"):
        self.parser = parser
        self.command_name = command_name
        self.help_text = help_text

        self.add_command(self.command_name, self.help_text)
    
    def add_command(self, command_name, help_text):
        self.subparser = self.parser.add_parser(command_name, help=help_text)

    def add_argument(self, arg, help_text = "Pas d'aide disponible", required = False, **kwargs):
        self.subparser.add_argument(arg, help=help_text, required=required, **kwargs)

    def match(self, command_name):
        return command_name == self.command_name

    def execute(self, namespace, **kwargs):
        pass
