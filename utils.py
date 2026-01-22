import requests
import os
import importlib.util
from cli_plugins.base.CLI_plugin import CLIPlugin

def get_mol_file(chebi_id):
    """
    Récupère le fichier mol d'un ID
    A implémenter : l'utilisation d'un cache et d'une BD
    
    L'endpoint utilisé est : /chebi/backend/api/public/molfile/{id}/
    
    :param chebi_id: L'ID Chebi de la molécule
    """

    url = f"https://www.ebi.ac.uk/chebi/backend/api/public/molfile/{chebi_id}"

    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return None
    
def get_all_plugins(path = "./cli_plugins"):
    # https://www.pythonmorsels.com/dynamically-importing-modules/
    # ^^^^^^^^^^^^ Pour le faire proprement, tah l'époque j'avais pas fait comme ça
    plugins = []
    for filename in os.listdir(path):
        if not filename.endswith(".py"): continue
        module_name = filename[:-3]
        module_path = os.path.join(path, filename)

        spec = importlib.util.spec_from_file_location(module_name, module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        for attr in dir(module):
            obj = getattr(module, attr)
            #on récup QUE la classe qui nous intéresse azy le reste
            if isinstance(obj, type) and issubclass(obj, module.CLIPlugin) and obj is not module.CLIPlugin:
                plugins.append(obj)
    return plugins

if __name__ == "__main__":
    #H20
    print(get_mol_file("15377")) 