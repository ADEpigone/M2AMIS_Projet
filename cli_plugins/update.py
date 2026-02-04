from cli_plugins.base.CLI_plugin import CLIPlugin
from utils import has_db_changed, get_mol_lite, prep_db_load, MOL_LITE_URL, MOL_LITE_UPDATE

class UpdatePlugin(CLIPlugin):

    def __init__(self, parser):
        super().__init__(parser, command_name="update", help_text = "Vérifie que la BD locale soit à jour")
        #self.add_argument("--bd_path", help_text="Path de la bd", required=True)

    def execute(self, namespace, **kwargs):
        print("Vérification de la mise à jour de la base de données locale...")
        info_date = {}
        if has_db_changed(MOL_LITE_URL, MOL_LITE_UPDATE, info_date):
            content = get_mol_lite(MOL_LITE_URL)
            prep_db_load(content, info_date)
        print("Base de données locale à jour.")