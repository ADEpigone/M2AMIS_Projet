import os
import csv
import random
import re
import importlib.util
import gzip
import io
import urllib.request

import requests
from rdkit import Chem, RDLogger
from tqdm import tqdm

from Chebi.ontology.parser import parse_obo
from Chebi.ontology.ontology_tree import OntologyTree
from graph import MoleculeGraph

MOL_LITE_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite.sdf.gz"
MOL_LITE_UPDATE = os.path.join("DB_updates", "mol_lite_last_update.txt")

ONTOLOGY_OBO_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"
ONTOLOGY_UPDATE = os.path.join("DB_updates", "ontology_last_update.txt")
ONTOLOGY_CACHE = os.path.join("DB_updates", "chebi_ontology.obo")

#ça devrait pas bougé c'est piqué de moleculenet
DEFAULT_ESOL_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv"


_ontology_tree_cache = None


def normalize_local_path(path: str) -> str:
    return os.path.normpath(str(path).replace("\\", os.sep))

def normalize_chebi_id(raw_value):
    value = str(raw_value).strip()
    if not value:
        return value

    upper = value.upper()
    if upper.startswith("CHEBI:"):
        suffix = value.split(":", 1)[1].strip()
        return f"CHEBI:{suffix}"
    if upper.startswith("CHEBI_"):
        suffix = value.split("_", 1)[1].strip()
        return f"CHEBI:{suffix}"

    chebi_tag_match = re.search(r"CHEBI[:_]\s*(\d+)", value, flags=re.IGNORECASE)
    if chebi_tag_match:
        return f"CHEBI:{chebi_tag_match.group(1)}"

    query_match = re.search(r"[?&]chebiId=([^&#]+)", value, flags=re.IGNORECASE)
    if query_match:
        return normalize_chebi_id(query_match.group(1))

    lower = value.lower()
    if "chebi" in lower:
        trailing_id_match = re.search(r"/(\d+)(?:[/?#].*)?$", value)
        if trailing_id_match:
            return f"CHEBI:{trailing_id_match.group(1)}"

    if value.isdigit():
        return f"CHEBI:{value}"

    return f"CHEBI:{value}"


def _write_remote_date(info_date):
    if not info_date:
        return
    if "local_path" not in info_date or "remote_date" not in info_date:
        return
    with open(info_date["local_path"], "w") as f:
        f.write(" ".join(info_date["remote_date"]))



def _select_log_s_key(fieldnames):
    candidates = [
        "logS",
        "log_s",
        "measured log solubility in mols per litre",
    ]
    for key in candidates:
        if key in fieldnames:
            return key
    return None


def download_esol(output_path=None, url=DEFAULT_ESOL_URL):
    if output_path is None:
        output_path = os.path.join("datasets", "esol.csv")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tmp_path = output_path + ".tmp"

    print(f"Telechargement ESOL depuis {url}")
    urllib.request.urlretrieve(url, tmp_path)

    with open(tmp_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise RuntimeError("CSV ESOL sans en-tete.")
        log_s_key = _select_log_s_key(reader.fieldnames)
        if log_s_key is None:
            raise RuntimeError("Colonne logS introuvable dans le CSV ESOL.")

        with open(output_path, "w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=["name", "smiles", "logS"])
            writer.writeheader()
            for row in reader:
                smiles = row.get("smiles") or row.get("SMILES")
                log_s = row.get(log_s_key)
                name = row.get("name") or row.get("Compound ID") or ""
                if smiles is None or log_s is None:
                    continue
                writer.writerow({"name": name, "smiles": smiles, "logS": log_s})

    os.remove(tmp_path)
    print(f"ESOL sauvegarde dans {output_path}")

def build_kernel(kernel_name: str, similarity: str):
    from similarites.builtin_similarity import BuiltinSimilarity
    from similarites.cwl_kernel import CWLKernel
    from similarites.ontology_similarity import OntologySimilarity
    if kernel_name == "cwl":
        return CWLKernel(similarity=similarity), True
    if kernel_name in {"morgan", "rdkit"}:
        return BuiltinSimilarity(fingerprint=kernel_name, similarity=similarity), True
    if kernel_name == "ontology":
        ontology = load_ontology()
        return OntologySimilarity(ontology), False
    raise ValueError(f"Kernel inconnu: {kernel_name}")

def load_esol_dataset(csv_path, limit=None, seed=42):
    if not os.path.exists(csv_path):
        print(f"Dataset ESOL introuvable: {csv_path}")
        return [], {"total": 0, "skipped": 0}

    molecules = []
    skipped = 0

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            if limit is not None and len(molecules) >= limit:
                break

            smiles = (row.get("smiles") or row.get("SMILES") or "").strip()
            log_s = row.get("logS")
            name = row.get("name") or row.get("Compound ID") or f"esol_{idx}"

            if not smiles or log_s is None:
                skipped += 1
                continue

            try:
                log_s_val = float(log_s)
                graph = MoleculeGraph.from_smiles(smiles, chebi_id=name)
                if graph is None or not graph.nodes:
                    raise ValueError("invalid graph")
            except Exception:
                skipped += 1
                continue

            molecules.append(
                {
                    "chebi_id": name,
                    "name": name,
                    "graph": graph,
                    "properties": {"logS": log_s_val, "smiles": smiles},
                }
            )

    random.Random(seed).shuffle(molecules)
    return molecules, {"total": len(molecules), "skipped": skipped, "source": "esol"}

def ensure_chebi_db_up_to_date():
    print("Vérification de la base locale ChEBI...")
    info_date = {}
    if has_db_changed(MOL_LITE_URL, MOL_LITE_UPDATE, info_date):
        content = get_mol_lite(MOL_LITE_URL)
        prep_db_load(content, info_date)
    print("Base locale ChEBI prête.")

def print_ids(id1, id2):
    chebi_id1 = normalize_chebi_id(id1)
    chebi_id2 = normalize_chebi_id(id2)
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

    normalized_local_path = normalize_local_path(local_db_date_path)

    response = requests.head(url)
    last_modified = response.headers.get("Last-Modified")
    if not last_modified:
        print("Header Last-Modified introuvable, update forcée.")
        return True
    remote_db_date = last_modified.split()[1:4]

    if info_date is not None:
        info_date['local_path'] = normalized_local_path
        info_date['remote_date'] = remote_db_date
    local_dir = os.path.dirname(normalized_local_path)
    if local_dir:
        os.makedirs(local_dir, exist_ok=True)

    if os.path.exists(normalized_local_path):
        with open(normalized_local_path, "r") as f:
            local_db_date = f.read().split()

        if remote_db_date == local_db_date:
            print("DB locale à jour !")
            return False
        print("Mise à jour de la DB locale requise !")
    else:
        print("Pas de DB locale !")

    if info_date is not None:
        info_date['local_path'] = normalized_local_path
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
                chebi_id = mol.GetProp("ChEBI ID") if mol.HasProp("ChEBI ID") else None
                name = mol.GetProp("ChEBI NAME") if mol.HasProp("ChEBI NAME") else None
                if name is None and mol.HasProp("_Name"):
                    name = mol.GetProp("_Name")
                mol_block = Chem.MolToMolBlock(mol)
            except Exception:
                continue

            batch.append((chebi_id, name, mol_block))
            
            if len(batch) >= 1000:
                db.update_table(batch)
                batch = []

        # Pour le dernier batch
        if batch:
            db.update_table(batch)
        pbar.close()

    _write_remote_date(info_date)
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

    normalized_cache_path = normalize_local_path(cache_path)
    cache_dir = os.path.dirname(normalized_cache_path)
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)

    with open(normalized_cache_path, 'w', encoding='utf-8') as f:
        f.write(obo_text)

    _write_remote_date(info_date)

    print("Ontologie ChEBI téléchargée et stockée !")
    return obo_text


def load_ontology(cache_path=ONTOLOGY_CACHE, force_download=False) -> OntologyTree:
    """
    Charge l'ontologie ChEBI.
    Télécharge automatiquement si pas en cache ou si mise à jour disponible.
    """
    normalized_cache_path = normalize_local_path(cache_path)

    global _ontology_tree_cache
    if _ontology_tree_cache is not None and not force_download:
        return _ontology_tree_cache

    obo_text = None

    # Vérifier si on doit mettre à jour
    info_date = {}
    try:
        needs_update = has_ontology_changed(info_date=info_date)
    except Exception:
        needs_update = not os.path.exists(normalized_cache_path)

    if needs_update or force_download or not os.path.exists(normalized_cache_path):
        obo_text = download_and_store_ontology(cache_path=normalized_cache_path, info_date=info_date)
    else:
        print("Chargement de l'ontologie depuis le cache local...")
        with open(normalized_cache_path, 'r', encoding='utf-8') as f:
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

    # Normaliser l'id (accepte "12345", "CHEBI:12345" ou URL ChEBI)
    chebi_id = normalize_chebi_id(chebi_id)

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
        module_name = f"_cli_plugin_{filename[:-3]}"
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

