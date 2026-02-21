import random

import numpy as np
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

from Chebi.CheBi2 import CheBi2
from analysis.clusters import DB_PATH
from graph import MoleculeGraph
from similarites.builtin_similarity import BuiltinSimilarity
from similarites.cwl_kernel import CWLKernel
from similarites.ontology_similarity import OntologySimilarity
from utils import load_ontology


TARGET_CLASSES = {
    "Amino acid": "CHEBI:33709",
    "Fatty acid": "CHEBI:35366",
    "Monosaccharide": "CHEBI:35381",
    "Steroid": "CHEBI:35341",
    "Flavonoid": "CHEBI:47916",
    "Alkaloid": "CHEBI:22315",
    "Porphyrin": "CHEBI:26214",
    "Lactone": "CHEBI:24973",
}


def compute_correlation_syntactic_semantic(molecules, syn_kernel, ontology_sim):
    print("\n" + "=" * 80)
    print("SYNTAXIQUE vs SÉMANTIQUE")
    print("=" * 80)

    n = len(molecules)
    if n < 2:
        print("Pas assez de molécules pour calculer une corrélation.")
        return

    print(f"\nCalcul des similarités pour {n} molécules ({n * (n - 1) // 2} paires)...")

    print("Calcul des fingerprints...")
    fingerprints = [syn_kernel.calculate_fingerprint(m["graph"]) for m in tqdm(molecules, desc="Fingerprints")]

    syntactic_sims = []
    semantic_sims = []

    print("Calcul des similarités...")
    for i in tqdm(range(n), desc="Avancement"):
        for j in range(i + 1, n):
            syntactic_sims.append(syn_kernel.calculate_similarity(fingerprints[i], fingerprints[j]))
            semantic_sims.append(ontology_sim.calculate_similarity(molecules[i]["graph"], molecules[j]["graph"]))

    syntactic_sims = np.array(syntactic_sims)
    semantic_sims = np.array(semantic_sims)

    pearson_r, pearson_p = pearsonr(syntactic_sims, semantic_sims)
    spearman_r, spearman_p = spearmanr(syntactic_sims, semantic_sims)

    print("\nRÉSULTATS:")
    print(f"  Corrélation de Pearson:   r = {pearson_r:.4f} (p = {pearson_p:.2e})")
    print(f"  Corrélation de Spearman:  ρ = {spearman_r:.4f} (p = {spearman_p:.2e})")

    print("\nSTATISTIQUES:")
    print(f"  Similarité syntaxique:   {syntactic_sims.mean():.3f} ± {syntactic_sims.std():.3f}")
    print(f"  Similarité sémantique:   {semantic_sims.mean():.3f} ± {semantic_sims.std():.3f}")

    print("\nCORRÉLATION PAR NIVEAUX DE SIMILARITÉ SYNTAXIQUE:")
    bins = [(0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
    for low, high in bins:
        mask = (syntactic_sims >= low) & (syntactic_sims < high)
        if mask.sum() > 0:
            avg_sem = semantic_sims[mask].mean()
            count = mask.sum()
            print(f"  Syntaxique [{low:.1f}-{high:.1f}]: {count:6d} paires / Sémantique moyen = {avg_sem:.3f}")

    print("=" * 80 + "\n")


def build_multiclass_sample(ontology, db, molecules_per_class=40, seed=42):
    random.seed(seed)
    selected = []

    print("Création d'un échantillon multi-classes...")
    for class_name, chebi_id in TARGET_CLASSES.items():
        node = ontology.get_node(chebi_id) if hasattr(ontology, "get_node") else None
        if not node:
            continue

        descendants = node.get_descendants() or set()
        descendants.add(node.chebi_id)
        candidates = list(descendants)
        random.shuffle(candidates)

        count = 0
        for cid in candidates:
            if count >= molecules_per_class:
                break

            mol_data = db.get_mol(cid)
            if not mol_data:
                continue

            try:
                graph = MoleculeGraph.from_moltext(mol_data, chebi_id=cid)
                if len(graph.nodes) <= 1:
                    continue
            except Exception:
                continue

            selected.append(
                {
                    "chebi_id": cid,
                    "name": f"{class_name}_{cid}",
                    "graph": graph,
                    "class_label": class_name,
                }
            )
            count += 1

        print(f"Ajouté {count} molécules de {class_name}.")

    return selected


def run_correlation_benchmark(molecules_per_class=40, kernel_similarity="cosine", fingerprint_method="cwl"):
    seed = 42
    db = CheBi2(DB_PATH)
    ontology = load_ontology()

    molecules = build_multiclass_sample(
        ontology=ontology,
        db=db,
        molecules_per_class=molecules_per_class,
        seed=seed,
    )

    print(f"\nLancement du calcul sur {len(molecules)} molécules variées...")
    if fingerprint_method == "cwl":
        syn_kernel = CWLKernel(similarity=kernel_similarity)
    elif fingerprint_method in {"morgan", "rdkit"}:
        syn_kernel = BuiltinSimilarity(fingerprint=fingerprint_method, similarity=kernel_similarity)
    else:
        raise ValueError("fingerprint_method doit être cwl, morgan ou rdkit")

    print(f"Fingerprint utilisé: {fingerprint_method}")
    ontology_sim = OntologySimilarity(ontology)
    compute_correlation_syntactic_semantic(molecules, syn_kernel, ontology_sim)


if __name__ == "__main__":
    run_correlation_benchmark()
