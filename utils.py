import requests
import os
import importlib.util
import gzip
import io
from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger

from Chebi.ontology.parser import parse_obo
from Chebi.ontology.ontology_tree import OntologyTree

MOL_LITE_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite.sdf.gz"
MOL_LITE_UPDATE = "DB_updates\\mol_lite_last_update.txt"

ONTOLOGY_OBO_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"
ONTOLOGY_UPDATE = "DB_updates\\ontology_last_update.txt"
ONTOLOGY_CACHE = "DB_updates\\chebi_ontology.obo"

_ontology_tree_cache = None

def print_ids(id1, id2):
    chebi_id1 = id1 if str(id1).upper().startswith("CHEBI:") else f"CHEBI:{id1}"
    chebi_id2 = id2 if str(id2).upper().startswith("CHEBI:") else f"CHEBI:{id2}"
    print(f"Molécule 1: https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id1}")
    print(f"Molécule 2: https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id2}")

def check_none(value, id):
    if value is None or value == "":
        print(f"L'ID {id} n'est pas une molécule valide pour CheBi")
        return True
    return False

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

    info_date['local_path'] = local_db_date_path
    info_date['remote_date'] = remote_db_date
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
    import Chebi.CheBi2 as Chebi2
    if response_content is None:
        raise Exception("Réponse vide")
    
    db = Chebi2.CheBi2("chebi2.db")

    RDLogger.DisableLog('rdApp.*')
    with gzip.open(io.BytesIO(response_content), 'rb') as f:
        # ForwardSDMolSupplier lit le flux compressé sans tout charger en RAM
        suppl = Chem.ForwardSDMolSupplier(f, sanitize=True, removeHs=False)
        
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


def has_ontology_changed(url=ONTOLOGY_OBO_URL, local_db_date_path=ONTOLOGY_UPDATE, info_date=None):
    """Vérifie si l'ontologie OBO distante a changé."""
    return has_db_changed(url=url, local_db_date_path=local_db_date_path, info_date=info_date)


def download_and_store_ontology(url=ONTOLOGY_OBO_URL, cache_path=ONTOLOGY_CACHE, info_date=None):
    """
    Télécharge le fichier chebi.obo.gz, le décompresse et le stocke localement.
    """
    print("Téléchargement de l'ontologie ChEBI...")
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Impossible de télécharger l'ontologie : HTTP {response.status_code}")

    with gzip.open(io.BytesIO(response.content), 'rb') as f:
        obo_text = f.read().decode('utf-8')

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, 'w', encoding='utf-8') as f:
        f.write(obo_text)

    # Mise à jour date locale
    if info_date and 'local_path' in info_date and 'remote_date' in info_date:
        with open(info_date['local_path'], "w") as f:
            f.write(" ".join(info_date['remote_date']))
        info_date = None

    print("Ontologie ChEBI téléchargée et stockée !")
    return obo_text


def load_ontology(cache_path=ONTOLOGY_CACHE, force_download=False) -> OntologyTree:
    """
    Charge l'ontologie ChEBI.
    Télécharge automatiquement si pas en cache ou si mise à jour disponible.
    """
    global _ontology_tree_cache
    if _ontology_tree_cache is not None and not force_download:
        return _ontology_tree_cache

    obo_text = None

    # Vérifier si on doit mettre à jour
    info_date = {}
    try:
        needs_update = has_ontology_changed(info_date=info_date)
    except Exception:
        needs_update = not os.path.exists(cache_path)

    if needs_update or force_download or not os.path.exists(cache_path):
        obo_text = download_and_store_ontology(info_date=info_date)
    else:
        print("Chargement de l'ontologie depuis le cache local...")
        with open(cache_path, 'r', encoding='utf-8') as f:
            obo_text = f.read()

    print("Parsing de l'ontologie...")
    _ontology_tree_cache = parse_obo(obo_text)
    print(f"Ontologie chargée : {_ontology_tree_cache}")
    return _ontology_tree_cache

def get_ontology_props(chebi_id: str, ontology: OntologyTree = None) -> dict:
    """
    Récupère les propriétés ontologiques d'un ID ChEBI.
    
    :return: dict avec name, definition, synonyms, parents, relationships, ancestors, etc.
    """
    if ontology is None:
        ontology = load_ontology()

    # Normaliser l'id (accepte "12345" ou "CHEBI:12345")
    if not chebi_id.startswith("CHEBI:"):
        chebi_id = f"CHEBI:{chebi_id}"

    node = ontology.get_node(chebi_id)
    if node is None:
        return {}

    return {
        "chebi_id": node.chebi_id,
        "name": node.name,
        "definition": node.definition,
        "synonyms": node.synonyms,
        "parents": [{"id": p.chebi_id, "name": p.name} for p in node.parents],
        "children": [{"id": c.chebi_id, "name": c.name} for c in node.children],
        "relationships": dict(node.properties),
        "ancestors": list(node.get_ancestors()),
        "depth": node.get_depth(),
    }

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

    ontology = load_ontology()
    
    props = get_ontology_props("16236", ontology)
    print(f"Nom: {props['name']}")
    print(f"Définition: {props['definition']}")
    print(f"Parents: {props['parents']}")
    print(f"Profondeur: {props['depth']}")
    print(f"Relationships: {props['relationships']}")
    print(f"Rôles : {props}")
    print(ontology.get_node('CHEBI:16236').get_roles())

    print(ontology.get_lca("CHEBI:15377", "CHEBI:16236"))