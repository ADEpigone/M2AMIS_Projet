import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

from analysis.clusters import DB_PATH, load_clusters
from utils import load_ontology
from Chebi.CheBi2 import CheBi2



def _safe_depth(node) -> int:
    """Retourne une profondeur si dispo, sinon 0 (fallback)."""
    if hasattr(node, "get_depth") and callable(getattr(node, "get_depth")):
        try:
            return int(node.get_depth())
        except Exception:
            return 0
    if hasattr(node, "depth"):
        try:
            return int(node.depth)
        except Exception:
            return 0
    return 0


def most_specific_family_id(ontology, mol_chebi_id: str, min_depth: int = 3) -> str | None:
    """
    Retourne UN identifiant de famille 'spécifique' pour la molécule.
    Heuristique:
      - candidats = {mol} U ancêtres(mol)
      - filtre profondeur >= min_depth (évite 'chemical entity', etc.)
      - prend le candidat avec profondeur maximale
    """
    mol_node = ontology.get_node(mol_chebi_id)
    if not mol_node:
        return None

    # candidats = mol + ancêtres
    cand_ids = set([mol_chebi_id])
    anc = mol_node.get_ancestors() or set()
    cand_ids |= set(anc)

    best_id = None
    best_depth = -1

    for cid in cand_ids:
        node = ontology.get_node(cid)
        if not node:
            continue
        d = _safe_depth(node)
        if d < min_depth:
            continue
        if d > best_depth:
            best_depth = d
            best_id = cid

    # Si rien ne passe le filtre, on fallback sur la molécule elle-même
    return best_id or mol_chebi_id


def is_ancestor_or_equal(ontology, ancestor_id: str, node_id: str) -> bool:
    """Vrai si ancestor_id est un ancêtre (ou égal) de node_id."""
    if ancestor_id == node_id:
        return True
    node = ontology.get_node(node_id)
    if not node:
        return False
    anc = node.get_ancestors() or set()
    return ancestor_id in anc


def dominant_family_ratio_for_cluster(
    ontology,
    cluster_members: list[dict],
    min_depth: int = 3
) -> tuple[float, str | None, Counter]:
    """
    Retourne:
      - ratio dominant (max_count / size)
      - id de la famille dominante (représentant)
      - distribution des familles dominantes candidates (Counter)
    """
    n = len(cluster_members)
    if n == 0:
        return 0.0, None, Counter()

    dominant_reps: list[str] = []   # liste des "familles dominantes candidates" (représentants)
    assign_counts = Counter()

    for m in cluster_members:
        mol_id = str(m.get("chebi_id"))
        fam_id = most_specific_family_id(ontology, mol_id, min_depth=min_depth)
        if fam_id is None:
            continue

        # Essayer d'assigner la molécule à un représentant existant
        assigned = False
        for rep in dominant_reps:
            # si la famille spécifique de la molécule est "parentée" (descendante)
            # d'une famille dominante déjà présente, on l'assigne à celle-ci.
            if is_ancestor_or_equal(ontology, rep, fam_id):
                assign_counts[rep] += 1
                assigned = True
                break

        if not assigned:
            dominant_reps.append(fam_id)
            assign_counts[fam_id] += 1

    if not assign_counts:
        return 0.0, None, Counter()

    dom_id, dom_count = assign_counts.most_common(1)[0]
    ratio = dom_count / n
    return ratio, dom_id, assign_counts


def dominant_ratios_for_all_clusters(
    ontology,
    cluster_map: dict,
    min_depth: int = 3,
    min_cluster_size: int = 2
):
    """
    Calcule DFR(C) pour tous les clusters en filtrant les petits clusters.
    """
    ratios = []
    details = {}

    for cid, members in cluster_map.items():
        if len(members) < min_cluster_size:
            continue

        ratio, dom_id, dist = dominant_family_ratio_for_cluster(
            ontology, members, min_depth=min_depth
        )
        ratios.append(ratio)
        details[cid] = {
            "size": len(members),
            "dominant_family_id": dom_id,
            "dominant_ratio": ratio,
            "distribution": dist,
        }

    return ratios, details


def plot_cumulative_curve(
    ratios: list[float],
    n_points: int = 101,
    title: str = "CDF des ratios de famille dominante",
    output_path: str = "cdf_dominant_family.png"
):
    """
    Sauvegarde y(t) = proportion de clusters avec ratio >= t
    dans un fichier image.
    """
    ratios = np.array(ratios, dtype=float)
    thresholds = np.linspace(0.0, 1.0, n_points)

    # Proportion de clusters qui dépassent le seuil
    y = [(ratios >= t).mean() for t in thresholds]

    plt.figure()
    plt.plot(thresholds, y)
    plt.xlabel("Seuil t sur le ratio dominant (DFR)")
    plt.ylabel("Proportion de clusters avec DFR ≥ t")
    plt.title(title)
    plt.grid(True)

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"[OK] Courbe cumulative sauvegardée dans : {output_path}")




if __name__ == "__main__":
    cluster_map = load_clusters("clusters_ontology_t0.7.json")
    ontology = load_ontology()

    chebi_db = CheBi2(DB_PATH)
    molecules = list(chebi_db.get_all_mols())

    # Calcule ratios par cluster
    ratios, details = dominant_ratios_for_all_clusters(ontology, cluster_map, min_depth=3)

    print(f"Nombre de clusters analysés: {len(ratios)}")
    print(f"Moyenne DFR: {np.mean(ratios):.3f} | Médiane DFR: {np.median(ratios):.3f}")
    print(f"% clusters avec DFR ≥ 0.7: {(np.mean(np.array(ratios) >= 0.7) * 100):.1f}%")

    # Trace la courbe cumulative
    plot_cumulative_curve(ratios, title="Courbe cumulative du ratio de famille dominante (clusters)", output_path="cdf_dominant_family_clusters.png")