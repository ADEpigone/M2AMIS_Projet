import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utils import load_ontology


def _safe_depth(node) -> int:
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


def load_clusters_from_json(json_path: str) -> dict[int, list[dict]]:
    path = Path(json_path)
    if not path.exists():
        raise FileNotFoundError(f"JSON introuvable: {json_path}")

    with open(path, "r") as f:
        data = json.load(f)

    clusters: dict[int, list[dict]] = {}
    for entry in data:
        cid = int(entry["cluster"])
        if cid not in clusters:
            clusters[cid] = []
        clusters[cid].append(entry)
    return clusters


def ancestors_with_self(ontology, mol_chebi_id: str) -> set[str]:
    node = ontology.get_node(mol_chebi_id)
    if not node:
        return set()
    anc = node.get_ancestors() or set()
    all_ids = set(anc)
    all_ids.add(mol_chebi_id)
    return all_ids


def dominant_consensus_for_cluster(
    ontology, cluster_members: list[dict]
) -> tuple[str | None, float, int, float, Counter]:
    """
    Dominance sur TOUT l'ensemble des ancetres:
      - on compte, pour chaque ancetre A, combien de molecules du cluster ont A
      - famille dominante = ancetre au meilleur score support-profondeur
      - ratio = support du dominant = count / n_molecules_valides
    Tie-break:
      - plus grand support
      - puis plus grande profondeur
      - puis id lexicographique (deterministe)
    """
    ancestor_counts = Counter()
    n_valid = 0

    for m in cluster_members:
        mol_id = str(m.get("chebi_id"))
        anc_ids = ancestors_with_self(ontology, mol_id)
        if not anc_ids:
            continue
        n_valid += 1
        for anc_id in anc_ids:
            ancestor_counts[anc_id] += 1

    if n_valid == 0 or not ancestor_counts:
        return None, 0.0, 0, 0.0, Counter()

    depths = {}
    max_depth = 0
    for anc_id in ancestor_counts.keys():
        node = ontology.get_node(anc_id)
        d = _safe_depth(node) if node else 0
        depths[anc_id] = d
        if d > max_depth:
            max_depth = d
    if max_depth <= 0:
        max_depth = 1

    best_id = None
    best_support = -1.0
    best_depth = -1
    best_select_score = -1.0
    for anc_id, count in ancestor_counts.items():
        depth = depths[anc_id]
        support = count / n_valid
        depth_norm = depth / max_depth
        select_score = support * depth_norm

        if select_score > best_select_score:
            best_id = anc_id
            best_support = support
            best_depth = depth
            best_select_score = select_score
        elif select_score == best_select_score:
            if support > best_support:
                best_id = anc_id
                best_support = support
                best_depth = depth
            elif support == best_support and depth > best_depth:
                best_id = anc_id
                best_depth = depth
            elif support == best_support and depth == best_depth and best_id is not None and anc_id < best_id:
                best_id = anc_id

    return best_id, best_support, best_depth, best_select_score, ancestor_counts


def consensus_scores_for_all_clusters(
    ontology,
    cluster_map: dict[int, list[dict]],
    min_cluster_size: int = 3,
    min_depth: int = 0,
) -> tuple[list[float], dict[int, dict]]:
    """
    Score consensus-profondeur (CDS):
      CDS = ratio_dominance * (profondeur_famille_dominante / profondeur_max_observee)
    """
    rows = []
    details = {}

    for cid, members in cluster_map.items():
        if len(members) < min_cluster_size:
            continue
        dom_id, ratio, depth, select_score, _ = dominant_consensus_for_cluster(ontology, members)
        if depth < min_depth:
            continue
        rows.append((cid, dom_id, ratio, depth, select_score, len(members)))

    if not rows:
        return [], details

    max_depth = max(r[3] for r in rows)
    if max_depth <= 0:
        max_depth = 1

    scores = []
    for cid, dom_id, ratio, depth, select_score, size in rows:
        depth_norm = depth / max_depth
        score = ratio * depth_norm
        scores.append(score)
        details[cid] = {
            "size": size,
            "dominant_family_id": dom_id,
            "dominance_ratio": ratio,
            "depth": depth,
            "depth_norm": depth_norm,
            "dominant_selection_score": select_score,
            "consensus_depth_score": score,
        }

    return scores, details


def plot_cumulative_curves(
    curves: dict[str, list[float]],
    n_points: int = 101,
    title: str = "CDF consensus-profondeur",
    output_path: str = "cdf_consensus_depth_comparison.png",
):
    thresholds = np.linspace(0.0, 1.0, n_points)

    plt.figure()
    for label, scores in curves.items():
        arr = np.array(scores, dtype=float)
        if len(arr) == 0:
            continue
        y = [(arr >= t).mean() for t in thresholds]
        plt.plot(thresholds, y, label=label)

    plt.xlabel("Seuil t sur CDS = ratio * profondeur_normalisee")
    plt.ylabel("Proportion de clusters avec CDS >= t")
    plt.title(title)
    plt.grid(True)
    plt.legend()

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[OK] Courbe sauvegardee: {output_path}")


def plot_cumulative_curve(
    scores: list[float],
    n_points: int = 101,
    title: str = "CDF consensus-profondeur",
    output_path: str = "cdf_consensus_depth.png",
):
    plot_cumulative_curves(
        curves={"CDS": scores},
        n_points=n_points,
        title=title,
        output_path=output_path,
    )


def plot_cds_vs_cluster_size(
    details: dict[int, dict],
    title: str,
    output_path: str,
):
    sizes = []
    scores = []
    for info in details.values():
        sizes.append(info["size"])
        scores.append(info["consensus_depth_score"])

    if not sizes:
        print(f"[WARN] Aucun point a tracer pour: {output_path}")
        return

    plt.figure()
    plt.scatter(sizes, scores, alpha=0.7, s=24)
    plt.xlabel("Taille du cluster")
    plt.ylabel("CDS")
    plt.title(title)
    plt.grid(True)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[OK] Nuage de points sauvegarde: {output_path}")


if __name__ == "__main__":
    ontology = load_ontology()

    cluster_configs = [
        {
            "label": "Morgan + Tanimoto",
            "json": "clusters_data_morgan_tanimoto.json",
            "tag": "morgan",
        },
        {
            "label": "CWL + Tanimoto",
            "json": "clusters_data_CWL_tanimoto.json",
            "tag": "cwl",
        },
    ]

    curves = {}
    for cfg in cluster_configs:
        label = cfg["label"]
        json_file = cfg["json"]
        tag = cfg["tag"]

        cluster_map = load_clusters_from_json(json_file)
        scores, details = consensus_scores_for_all_clusters(ontology, cluster_map, min_cluster_size=2)

        if len(scores) == 0:
            print(f"{label} | aucun cluster analyse.")
            continue

        curves[label] = scores
        arr = np.array(scores, dtype=float)
        print(f"{label} | Clusters analyses: {len(arr)}")
        print(f"{label} | Moyenne CDS: {arr.mean():.3f} | Mediane CDS: {np.median(arr):.3f}")
        print(f"{label} | % clusters avec CDS >= 0.5: {(np.mean(arr >= 0.5) * 100):.1f}%")

        plot_cds_vs_cluster_size(
            details,
            title=f"CDS en fonction de la taille des clusters ({label})",
            output_path=f"scatter_cds_vs_size_{tag}.png",
        )

    if not curves:
        raise RuntimeError("Aucune courbe a tracer. Verifie les JSON d'entree.")

    plot_cumulative_curves(
        curves,
        title="Comparaison des courbes cumulatives du consensus-profondeur",
    )
