import csv
import json
import os
import random
from datetime import datetime
from math import sqrt
import time
import numpy as np

from similarites.cwl_kernel import CWLKernel
from graph import MoleculeGraph
from similarites.builtin_similarity import BuiltinSimilarity

def load_esol_dataset(csv_path, limit=None, seed=42):
    if not os.path.exists(csv_path):
        print(f"Dataset ESOL introuvable: {csv_path}")
        return [], {"total": 0, "skipped": 0}

    rng = random.Random(seed)
    molecules = []
    skipped = 0

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            if limit is not None and len(molecules) >= limit:
                break
            smiles = row.get("smiles") or row.get("SMILES")
            if smiles is not None:
                smiles = smiles.strip()
            log_s = row.get("logS")
            name = row.get("name") or row.get("Compound ID") or f"esol_{idx}"
            if smiles is None or log_s is None:
                skipped += 1
                continue
            try:
                log_s_val = float(log_s)
            except (TypeError, ValueError):
                skipped += 1
                continue

            graph = MoleculeGraph.from_smiles(smiles, chebi_id=name)
            if graph is None or not graph.nodes:
                skipped += 1
                continue

            molecules.append({
                "chebi_id": name,
                "name": name,
                "graph": graph,
                "properties": {
                    "logS": log_s_val,
                    "smiles": smiles,
                },
            })

    rng.shuffle(molecules)
    info = {
        "total": len(molecules),
        "skipped": skipped,
        "source": "esol",
    }
    return molecules, info


def _sample_indices(n, sample_frac, seed):
    if sample_frac >= 1.0:
        return list(range(n))
    k = max(10, int(round(n * sample_frac)))
    k = min(k, n)
    rng = random.Random(seed)
    return rng.sample(range(n), k)


def _random_neighbor_indices(n, seed):
    rng = random.Random(seed)
    indices = []
    for i in range(n):
        choices = [x for x in range(n) if x != i]
        indices.append(rng.choice(choices))
    return indices


def _build_trials(n, sample_seeds, rand_seeds, sample_frac):
    trials = []
    for sseed in sample_seeds:
        indices = _sample_indices(n, sample_frac, sseed)
        for rseed in rand_seeds:
            rand_idx = _random_neighbor_indices(len(indices), rseed)
            trials.append({
                "indices": indices,
                "rand_idx": rand_idx,
                "sample_seed": sseed,
                "rand_seed": rseed,
            })
    return trials


def _mean_local_error_subset(valid_mols, vals, kernel, indices, rand_idx):
    subset_mols = [valid_mols[i] for i in indices]
    subset_vals = [vals[i] for i in indices]
    n = len(subset_mols)

    fingerprints = [kernel.calculate_fingerprint(m["graph"]) for m in subset_mols]
    deltas_nearest = []
    deltas_random = []

    # Precompute similarity matrix using symmetry to halve computations.
    sims = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            sim = kernel.calculate_similarity(fingerprints[i], fingerprints[j])
            sims[i][j] = sim
            sims[j][i] = sim

    for i in range(n):
        best_sim = -1
        best_idx = -1
        row = sims[i]
        for j in range(n):
            if i == j:
                continue
            sim = row[j]
            if sim > best_sim:
                best_sim = sim
                best_idx = j
        if best_idx != -1:
            deltas_nearest.append(abs(subset_vals[i] - subset_vals[best_idx]))
        deltas_random.append(abs(subset_vals[i] - subset_vals[rand_idx[i]]))

    mean_struct = float(np.mean(deltas_nearest)) if deltas_nearest else float("inf")
    mean_rand = float(np.mean(deltas_random)) if deltas_random else float("inf")
    return mean_struct, mean_rand


def _compute_stats(values):
    arr = np.array(values, dtype=float)
    n = len(arr)
    mean = float(arr.mean()) if n > 0 else float("nan")
    return {
        "n": n,
        "mean": mean,
    }


def _evaluate_kernel(valid_mols, vals, kernel, trials):
    struct_errors = []
    rand_errors = []
    ratios = []
    for trial in trials:
        mean_struct, mean_rand = _mean_local_error_subset(
            valid_mols,
            vals,
            kernel,
            trial["indices"],
            trial["rand_idx"],
        )
        struct_errors.append(mean_struct)
        rand_errors.append(mean_rand)
        ratio = mean_rand / mean_struct if mean_struct > 0 else 0.0
        ratios.append(ratio)

    return {
        "struct": _compute_stats(struct_errors),
        "rand": _compute_stats(rand_errors),
        "ratio": _compute_stats(ratios),
    }


def compare_similarity_kernels(
    molecules,
    prop,
    experiment_name,
    similarity_names,
    builtin_fingerprints,
    sample_frac=0.8,
    output_dir="clustering/results",
):
    valid_mols = []
    vals = []
    for m in molecules:
        props = m.get("properties")
        try:
            val = float(props.get(prop))
        except (TypeError, ValueError):
            continue
        vals.append(val)
        valid_mols.append(m)
    sample_seeds = [11, 23, 37, 51, 73]
    rand_seeds = [101, 131, 151, 173, 199]

    trials = _build_trials(len(valid_mols), sample_seeds, rand_seeds, sample_frac)

    results = {
        "experiment_name": experiment_name,
        "property": prop,
        "num_molecules": len(valid_mols),
        "num_trials": len(trials),
        "sample_frac": sample_frac,
        "sample_seeds": sample_seeds,
        "rand_seeds": rand_seeds,

        "comparisons": [],
    }

    for sim_name in similarity_names:
        print("\n" + "=" * 70)
        print(f"SIM={sim_name} | PROP={prop}")
        print("=" * 70)

        builtin_stats_by_fp = {}
        random_stats = None
        for fp_name in builtin_fingerprints:
            print(f"BUILTIN={fp_name}")
            builtin = BuiltinSimilarity(fingerprint=fp_name, similarity=sim_name)
            t1 = time.time()
            builtin_stats = _evaluate_kernel(valid_mols, vals, builtin, trials)
            t2 = time.time()
            builtin_stats_by_fp[fp_name] = builtin_stats
            if random_stats is None:
                random_stats = builtin_stats["rand"]
            print(f"Builtin kernel evaluation in {t2 - t1:.2f}s.")
            print(
                f"Baseline {fp_name} mean err: {builtin_stats['struct']['mean']:.3f} "
                f"(ratio {builtin_stats['ratio']['mean']:.2f})"
            )

        if random_stats is not None:
            print(
                f"Random baseline mean err: {random_stats['mean']:.3f}"
            )

        kernel = CWLKernel(similarity=sim_name)
        t1 = time.time()
        cwl_stats = _evaluate_kernel(valid_mols, vals, kernel, trials)
        t2 = time.time()
        print(f"CWL evaluation in {t2 - t1:.2f}s.")
        print(
            f"CWL -> mean err {cwl_stats['struct']['mean']:.3f} "
            f"(ratio {cwl_stats['ratio']['mean']:.2f})"
        )

        for fp_name, builtin_stats in builtin_stats_by_fp.items():
            results["comparisons"].append({
                "similarity": sim_name,
                "builtin_fingerprint": fp_name,
                "builtin_stats": builtin_stats,
                "cwl_stats": cwl_stats,
                "delta_vs_builtin": cwl_stats["struct"]["mean"] - builtin_stats["struct"]["mean"],
            })

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(output_dir, f"{experiment_name}_tuning_{timestamp}.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=True, indent=2)

    print(f"\nResultats sauvegardes dans {out_path}")
    return out_path


if __name__ == "__main__":
    esol_csv = os.path.join("datasets", "esol.csv")
    print("CHARGEMENT DATASET ESOL (logS)")
    molecules, dataset_info = load_esol_dataset(esol_csv, limit=None, seed=42)
    print(f"Dataset ESOL: {dataset_info}")
    if not molecules:
        raise SystemExit("Dataset ESOL vide. Lance datasets/download_esol.py pour le telecharger.")

    similarity_names = ["tanimoto", "dice", "cosine"]
    builtin_fingerprints = ["morgan", "rdkit"]

    compare_similarity_kernels(
        molecules,
        prop="logS",
        experiment_name="esol_logS",
        similarity_names=similarity_names,
        builtin_fingerprints=builtin_fingerprints,
        sample_frac=0.8,
    )
