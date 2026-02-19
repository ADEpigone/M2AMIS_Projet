import csv
import json
import os
import random
from datetime import datetime
import time
import numpy as np
from rdkit import DataStructs

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


def _fp_to_numpy(fp):
    arr = np.zeros((fp.GetNumBits(),), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(np.float32)


def _build_xy_for_fingerprint(molecules, fp_name):
    if fp_name == "cwl":
        kernel = CWLKernel(similarity="tanimoto")
    else:
        kernel = BuiltinSimilarity(fingerprint=fp_name, similarity="tanimoto")

    X = []
    y = []
    skipped = 0

    for m in molecules:
        props = m.get("properties", {})
        try:
            target = float(props.get("logS"))
        except (TypeError, ValueError):
            skipped += 1
            continue

        try:
            fp = kernel.calculate_fingerprint(m["graph"])
            vec = _fp_to_numpy(fp)
        except Exception:
            skipped += 1
            continue

        X.append(vec)
        y.append(target)

    if not X:
        return np.empty((0, 0), dtype=np.float32), np.empty((0,), dtype=np.float32), skipped

    return np.vstack(X), np.array(y, dtype=np.float32), skipped


def _agg_metric(values):
    arr = np.array(values, dtype=float)
    return {
        "n": int(len(arr)),
        "mean": float(arr.mean()) if len(arr) > 0 else float("nan"),
        "std": float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
    }


def _evaluate_regressor_cv(X, y, splits, model_name):
    from sklearn.dummy import DummyRegressor
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    if model_name == "dummy":
        build = lambda: DummyRegressor(strategy="mean")
    elif model_name == "ridge":
        build = lambda: Ridge(alpha=1.0)
    elif model_name == "rf":
        build = lambda: RandomForestRegressor(
            n_estimators=500,
            random_state=42,
            n_jobs=-1,
        )
    else:
        raise ValueError(f"Modele inconnu: {model_name}")

    folds = []
    mae_vals = []
    rmse_vals = []
    r2_vals = []

    for fold_id, (tr, te) in enumerate(splits, start=1):
        model = build()
        model.fit(X[tr], y[tr])
        pred = model.predict(X[te])

        mae = float(mean_absolute_error(y[te], pred))
        rmse = float(np.sqrt(mean_squared_error(y[te], pred)))
        r2 = float(r2_score(y[te], pred))

        folds.append({
            "fold": fold_id,
            "mae": mae,
            "rmse": rmse,
            "r2": r2,
        })
        mae_vals.append(mae)
        rmse_vals.append(rmse)
        r2_vals.append(r2)

    return {
        "folds": folds,
        "mae": _agg_metric(mae_vals),
        "rmse": _agg_metric(rmse_vals),
        "r2": _agg_metric(r2_vals),
    }


def compare_prediction_benchmark(
    molecules,
    experiment_name="esol_prediction",
    fingerprint_names=("morgan", "rdkit", "cwl"),
    model_names=("dummy", "ridge", "rf"),
    n_splits=5,
    output_dir="clustering/results",
):
    from sklearn.model_selection import KFold

    results = {
        "experiment_name": experiment_name,
        "target": "logS",
        "n_folds": n_splits,
        "fingerprint_results": [],
    }

    for fp_name in fingerprint_names:
        print("\n" + "=" * 70)
        print(f"PREDICTION | FP={fp_name} | TARGET=logS")
        print("=" * 70)

        t0 = time.time()
        X, y, skipped = _build_xy_for_fingerprint(molecules, fp_name)
        if len(X) == 0:
            print(f"Aucun echantillon valide pour {fp_name}.")
            continue

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        splits = list(kf.split(X, y))

        fp_result = {
            "fingerprint": fp_name,
            "n_samples": int(len(y)),
            "n_features": int(X.shape[1]),
            "skipped": int(skipped),
            "models": {},
        }

        for model_name in model_names:
            model_stats = _evaluate_regressor_cv(X, y, splits, model_name)
            fp_result["models"][model_name] = model_stats
            print(
                f"{model_name:>5} | "
                f"MAE={model_stats['mae']['mean']:.3f}±{model_stats['mae']['std']:.3f} | "
                f"RMSE={model_stats['rmse']['mean']:.3f}±{model_stats['rmse']['std']:.3f} | "
                f"R2={model_stats['r2']['mean']:.3f}±{model_stats['r2']['std']:.3f}"
            )

        fp_result["elapsed_sec"] = float(time.time() - t0)
        results["fingerprint_results"].append(fp_result)

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(output_dir, f"{experiment_name}_tuning_{timestamp}.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=True, indent=2)

    print(f"\nResultats prediction sauvegardes dans {out_path}")
    return out_path


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

    compare_prediction_benchmark(
        molecules,
        experiment_name="esol_prediction",
        fingerprint_names=("morgan", "rdkit", "cwl"),
        model_names=("dummy", "ridge", "rf"),
        n_splits=5,
    )
