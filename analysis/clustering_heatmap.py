import random

import matplotlib.pyplot as plt
import numpy as np
from rdkit import RDLogger
from scipy.cluster.hierarchy import linkage, leaves_list
from tqdm import tqdm

from Chebi.CheBi2 import CheBi2
from graph import MoleculeGraph

RDLogger.DisableLog("rdApp.*")

DB_PATH = "chebi2.db"


def _load_sampled_molecules(max_molecules: int | None = None) -> list[dict]:
    db = CheBi2(DB_PATH)

    molecules = []
    raw_data = [m for m in db.get_all_mols() if m[2]]
    if max_molecules is not None and max_molecules > 0:
        sample_size = min(int(max_molecules), len(raw_data))
        raw_data = random.sample(raw_data, sample_size)

    for chebi_id, mol_name, mol_file in tqdm(raw_data, total=len(raw_data), desc="Parsing"):
        if not mol_file:
            continue
        try:
            graph = MoleculeGraph.from_moltext(mol_file, chebi_id=chebi_id)
            if len(graph.nodes) == 0:
                continue
            molecules.append(
                {
                    "chebi_id": chebi_id,
                    "name": mol_name,
                    "graph": graph,
                }
            )
        except Exception:
            continue

        if max_molecules is not None and max_molecules > 0 and len(molecules) >= int(max_molecules):
            break

    return molecules


def generate_hierarchical_similarity_heatmap(
    sim_kernel,
    has_fingerprint: bool,
    dist_threshold: float,
    max_molecules: int | None,
    output_path: str = "heatmap_hierarchical_similarity.png",
):
    molecules = _load_sampled_molecules(max_molecules=max_molecules)
    n = len(molecules)

    if n < 2:
        raise RuntimeError("Pas assez de molécules.")

    print(f"{n} molécules prêtes pour la heatmap.")

    print("Calcul des fingerprints...")
    fingerprints = None
    if has_fingerprint:
        fingerprints = [sim_kernel.calculate_fingerprint(m["graph"]) for m in tqdm(molecules)]

    print("Calcul matrice de similarité / distance...")
    sim_matrix = np.eye(n, dtype=float)
    condensed_dist = []

    for i in tqdm(range(n), desc="Similarités"):
        for j in range(i + 1, n):
            if has_fingerprint:
                s = sim_kernel.calculate_similarity(fingerprints[i], fingerprints[j])
            else:
                s = sim_kernel.calculate_similarity(molecules[i]["graph"], molecules[j]["graph"])
            s = float(max(0.0, min(1.0, s)))
            sim_matrix[i, j] = s
            sim_matrix[j, i] = s
            condensed_dist.append(max(0.0, 1.0 - s))

    dist_array = np.array(condensed_dist, dtype=float)
    Z = linkage(dist_array, method="average")
    ordered_idx = leaves_list(Z)

    ordered_sim = sim_matrix[np.ix_(ordered_idx, ordered_idx)]

    plt.figure(figsize=(8, 7))
    im = plt.imshow(ordered_sim, cmap="viridis", vmin=0.0, vmax=1.0, interpolation="nearest")
    plt.colorbar(im, label="Similarité")
    plt.title(f"Heatmap des similarités (seuil={dist_threshold:.2f})")

    if n <= 60:
        labels = [str(molecules[i]["chebi_id"]) for i in ordered_idx]
        plt.xticks(range(n), labels, rotation=90, fontsize=6)
        plt.yticks(range(n), labels, fontsize=6)
    else:
        plt.xticks([])
        plt.yticks([])

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Heatmap sauvegardée: {output_path}")
