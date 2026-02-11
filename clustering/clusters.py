import numpy as np
import json
import random
from tqdm import tqdm
from pathlib import Path
import sys
sys.path.append("./..")
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from rdkit import RDLogger

from Chebi.CheBi2 import CheBi2
from graph import MoleculeGraph
from datasets.base.dataset import MoleculeEntry
from similarites.cwl_kernel import CWLKernel
# from similarites.builtin_similarity import BuiltinSimilarity
from utils import *
RDLogger.DisableLog('rdApp.*')

DB_PATH = "chebi2.db"
JSON_OUTPUT = "clusters_data.json"
MAX_MOLECULES = 10000
DIST_THRESHOLD = 0.75 # Seuil pour le clustering

def run_clustering_and_save():
    db = CheBi2(DB_PATH)

    molecules = []
    raw_data = [m for m in db.get_all_mols() if m[2]]
    raw_data = random.sample(raw_data, MAX_MOLECULES)
    for chebi_id, mol_name, mol_file in tqdm(raw_data, total=len(raw_data), desc="Parsing"):
        if not mol_file: continue
        try:
            graph = MoleculeGraph.from_moltext(mol_file, chebi_id=chebi_id)
            if len(graph.nodes) == 0: continue

            molecules.append({
                "chebi_id": chebi_id,
                "name": mol_name,
                "graph": graph, 
            })
        except Exception:
            continue

        if MAX_MOLECULES and len(molecules) >= MAX_MOLECULES:
            break
    
    print(f"{len(molecules)} molécules prêtes.")

    sim_kernel = CWLKernel(similarity="tanimoto")
    
    print("Calcul des fingerprints...")
    fingerprints = [sim_kernel.calculate_fingerprint(m["graph"]) for m in tqdm(molecules)]

    print("Calcul matrice de distance...")
    n = len(molecules)
    condensed_dist = []
    
    for i in tqdm(range(n), desc="Distances"):
        for j in range(i + 1, n):
            s = sim_kernel.calculate_similarity(fingerprints[i], fingerprints[j])
            condensed_dist.append(max(0.0, 1.0 - s))
    
    dist_array = np.array(condensed_dist)

    print("Clustering...")
    Z = linkage(dist_array, method='average')
    labels = fcluster(Z, t=DIST_THRESHOLD, criterion='distance')

    export_data = []
    for i, m in enumerate(molecules):
        export_data.append({
            "chebi_id": m["chebi_id"],
            "name": m["name"],
            "cluster": int(labels[i])
        })

    with open(JSON_OUTPUT, "w") as f:
        json.dump(export_data, f, indent=4)
    
    print(f"Fichier save dans : {JSON_OUTPUT}")


def load_clusters():
    """
    Charge juste le JSON pour analyse rapide sans tout recalculer.
    """
    info_date = {}
    if has_db_changed(MOL_LITE_URL, MOL_LITE_UPDATE, info_date):
        content = get_mol_lite(MOL_LITE_URL)
        prep_db_load(content, info_date)

    chebi_db = CheBi2(DB_PATH)

    if not Path(JSON_OUTPUT).exists():
        print("Fichier JSON introuvable.")
        return

    with open(JSON_OUTPUT, "r") as f:
        data = json.load(f)

    for elem in data:
        elem["graph"] = MoleculeGraph.from_moltext(chebi_db.get_mol(elem["chebi_id"]))
    
    clusters = {}
    for entry in data:
        c_id = entry["cluster"]
        if c_id not in clusters: clusters[c_id] = []
        clusters[c_id].append(entry)

    print(f"Nombre de clusters trouvés : {len(clusters)}")

    return clusters


if __name__ == "__main__":

    run_clustering_and_save()
    
    #load_clusters()