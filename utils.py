import requests
import os
import importlib.util
import gzip
import io
import Chebi.CheBi2 as Chebi2
from tqdm import tqdm
from cli_plugins.base.CLI_plugin import CLIPlugin
from rdkit import Chem
from rdkit import RDLogger

MOL_LITE_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite.sdf.gz"
MOL_LITE_UPDATE = "DB_updates\\mol_lite_last_update.txt"

def get_mol_lite(url=MOL_LITE_URL):
    """
    Télécharge le fichier chebi_lite.sdf.gz

    contenu :
    - mol_file
    - id
    - name
    - star_lvl

    :return: Le contenu du fichier SDF compressé
    """
    print("Téléchargement du fichier")
    response = requests.get(url)
    if response.status_code == 200:
        return response.content
    else:
        return None

def has_db_changed(url=None, local_db_date_path=None, info_date=None):
    if url is None or local_db_date_path is None:
        raise Exception(f"URL et/ou DB_file_version manquant ({url=}, {local_db_date_path=})")
    response = requests.head(url)
    remote_db_date = response.headers.get('Last-Modified').split()[1:4]

    os.makedirs(os.path.dirname(local_db_date_path), exist_ok=True)
    if os.path.exists(local_db_date_path):
        with open(local_db_date_path, "r") as f:
            local_db_date = f.read().split()

        if remote_db_date == local_db_date:
            print("DB locale à jour !")
            return False
        print("Mise à jour de la DB locale requise !")
    else:
        print("Pas de DB locale !")

    info_date['local_path'] = local_db_date_path
    info_date['remote_date'] = remote_db_date
    return True

def prep_db_load(response_content, info_date=None):
    if response_content is None:
        raise Exception("Réponse vide")
    
    db = Chebi2.CheBi2("chebi2.db")

    RDLogger.DisableLog('rdApp.*')
    with gzip.open(io.BytesIO(response_content), 'rb') as f:
        # ForwardSDMolSupplier lit le flux compressé sans tout charger en RAM
        suppl = Chem.ForwardSDMolSupplier(f, sanitize=False, removeHs=False)
        
        batch = []
        pbar = tqdm(desc="Parsing ChEBI_SDF_file", unit=" mol")
        for mol in suppl:
            pbar.update(1)
            if mol is None: continue  # skip les mol mal formées
            
            try:
                chebi_id = mol.GetProp('ChEBI ID') if mol.HasProp('ChEBI ID') else None
                name = mol.GetProp('ChEBI NAME') if mol.HasProp('ChEBI NAME') else None
                mol_format = Chem.MolToMolBlock(mol) if mol else None
            except Exception as e:
                continue
            batch.append((chebi_id, name, mol_format))
            
            if len(batch) >= 1000:
                db.update_table(batch)
                batch = []

        # Pour le dernier batch
        if batch:
            db.update_table(batch)
        pbar.close()

    # Mise à jour date locale
    with open(info_date['local_path'], "w") as f:
        f.write(" ".join(info_date['remote_date']))
    info_date = None
    print("DB Chebi2 Up To Date !")

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
    #print(get_mol_file("15377")) 
    #print(has_db_changed(MOL_LITE_URL, MOL_LITE_UPDATE))
    info_date = {}
    if has_db_changed(MOL_LITE_URL, MOL_LITE_UPDATE, info_date):
        content = get_mol_lite(MOL_LITE_URL)
        prep_db_load(content, info_date)