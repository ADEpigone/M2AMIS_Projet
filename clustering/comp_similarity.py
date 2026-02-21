import random

import numpy as np

from similarites.builtin_similarity import BuiltinSimilarity
from similarites.cwl_kernel import CWLKernel


def _sample_indices(n, sample_frac, seed):
    if sample_frac >= 1.0:
        return list(range(n))
    k = max(10, int(round(n * sample_frac)))
    k = min(k, n)
    return random.Random(seed).sample(range(n), k)


def _random_neighbor_indices(n, seed):
    if n <= 1:
        return [0] * max(0, n)
    rng = random.Random(seed)
    return [rng.choice([x for x in range(n) if x != i]) for i in range(n)]


def _mean_local_error_subset(valid_mols, values, kernel, indices, random_indices):
    if len(indices) <= 1:
        return float("inf"), float("inf"), 0.0

    subset_mols = [valid_mols[i] for i in indices]
    subset_vals = [values[i] for i in indices]
    n = len(subset_mols)

    fps = [kernel.calculate_fingerprint(m["graph"]) for m in subset_mols]

    deltas_nearest = []
    deltas_random = []

    for i in range(n):
        best_sim = -1.0
        best_idx = -1
        for j in range(n):
            if i == j:
                continue
            sim = kernel.calculate_similarity(fps[i], fps[j])
            if sim > best_sim:
                best_sim = sim
                best_idx = j

        if best_idx != -1:
            deltas_nearest.append(abs(subset_vals[i] - subset_vals[best_idx]))
        deltas_random.append(abs(subset_vals[i] - subset_vals[random_indices[i]]))

    mean_struct = float(np.mean(deltas_nearest)) if deltas_nearest else float("inf")
    mean_rand = float(np.mean(deltas_random)) if deltas_random else float("inf")
    ratio = mean_rand / mean_struct if mean_struct > 0 else 0.0

    return mean_struct, mean_rand, ratio


def _evaluate_kernel(valid_mols, values, kernel, sample_frac, sample_seeds, rand_seeds):
    struct_vals = []
    rand_vals = []
    ratio_vals = []

    for sample_seed in sample_seeds:
        indices = _sample_indices(len(valid_mols), sample_frac, sample_seed)
        for rand_seed in rand_seeds:
            random_indices = _random_neighbor_indices(len(indices), rand_seed)
            s, r, ratio = _mean_local_error_subset(valid_mols, values, kernel, indices, random_indices)
            struct_vals.append(s)
            rand_vals.append(r)
            ratio_vals.append(ratio)

    return (
        float(np.mean(struct_vals)) if struct_vals else float("inf"),
        float(np.mean(rand_vals)) if rand_vals else float("inf"),
        float(np.mean(ratio_vals)) if ratio_vals else 0.0,
    )


def compare_similarity_kernels(
    molecules,
    prop,
    experiment_name,
    similarity_names,
    builtin_fingerprints,
    sample_frac=0.8,
):
    print(f"\nLancement benchmark similarité: {experiment_name}")
    valid_mols = []
    values = []
    for molecule in molecules:
        try:
            value = float(molecule["properties"][prop])
        except Exception:
            continue
        valid_mols.append(molecule)
        values.append(value)

    if len(valid_mols) < 2:
        print("Pas assez de molécules valides pour benchmark similarité (minimum: 2).")
        return None

    sample_seeds = [11, 23, 37]
    rand_seeds = [101, 131, 151]

    for sim_name in similarity_names:
        print("\n" + "=" * 70)
        print(f"SIM={sim_name} | PROP={prop}")
        print("=" * 70)

        cwl = CWLKernel(similarity=sim_name)
        cwl_mean_struct_err, cwl_mean_rand_err, cwl_mean_ratio = _evaluate_kernel(
            valid_mols,
            values,
            cwl,
            sample_frac,
            sample_seeds,
            rand_seeds,
        )

        print(
            f"cwl | mean={cwl_mean_struct_err:.3f} "
            f"(rand={cwl_mean_rand_err:.3f}, ratio={cwl_mean_ratio:.2f})"
        )

        for fp_name in builtin_fingerprints:
            try:
                builtin = BuiltinSimilarity(fingerprint=fp_name, similarity=sim_name)
            except ValueError:
                print(f"{fp_name} ignoré (fingerprint non supporté pour builtin).")
                continue

            builtin_mean_struct_err, builtin_mean_rand_err, builtin_mean_ratio = _evaluate_kernel(
                valid_mols,
                values,
                builtin,
                sample_frac,
                sample_seeds,
                rand_seeds,
            )

            print(
                f"{fp_name} | "
                f"mean={builtin_mean_struct_err:.3f} "
                f"(rand={builtin_mean_rand_err:.3f}, ratio={builtin_mean_ratio:.2f})"
            )

    print("\nBenchmark similarité terminé.")
    return None
