import numpy as np
from rdkit import DataStructs

from similarites.builtin_similarity import BuiltinSimilarity
from similarites.cwl_kernel import CWLKernel
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold


def _fp_to_numpy(fp):
    arr = np.zeros((fp.GetNumBits(),), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(np.float32)


def _build_xy_for_fingerprint(molecules, fp_name):
    kernel = CWLKernel(similarity="tanimoto") if fp_name == "cwl" else BuiltinSimilarity(fingerprint=fp_name, similarity="tanimoto")

    x_values = []
    y_values = []
    skipped = 0

    for molecule in molecules:
        try:
            target = float(molecule["properties"]["logS"])
            fp = kernel.calculate_fingerprint(molecule["graph"])
            vec = _fp_to_numpy(fp)
        except Exception:
            skipped += 1
            continue

        x_values.append(vec)
        y_values.append(target)

    if not x_values:
        return np.empty((0, 0), dtype=np.float32), np.empty((0,), dtype=np.float32), skipped

    return np.vstack(x_values), np.array(y_values, dtype=np.float32), skipped


def _evaluate_model_cv(x_data, y_data, n_splits, model_name):

    if model_name == "ridge":
        model_factory = lambda: Ridge(alpha=1.0)
    elif model_name == "rf":
        model_factory = lambda: RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
    else:
        raise ValueError(f"Modele inconnu: {model_name}")

    mae_values = []
    r2_values = []

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    for train_idx, test_idx in kf.split(x_data, y_data):
        model = model_factory()
        model.fit(x_data[train_idx], y_data[train_idx])
        pred = model.predict(x_data[test_idx])

        mae_values.append(float(mean_absolute_error(y_data[test_idx], pred)))
        r2_values.append(float(r2_score(y_data[test_idx], pred)))

    return {
        "mae_mean": float(np.mean(mae_values)),
        "mae_std": float(np.std(mae_values, ddof=1)) if len(mae_values) > 1 else 0.0,
        "r2_mean": float(np.mean(r2_values)),
        "r2_std": float(np.std(r2_values, ddof=1)) if len(r2_values) > 1 else 0.0,
    }


def compare_prediction_benchmark(
    molecules,
    experiment_name="esol_prediction",
    fingerprint_names=("morgan", "rdkit", "cwl"),
    model_names=("ridge", "rf"),
    n_splits=5
):
    print(f"\nLancement benchmark: {experiment_name}")

    for fp_name in fingerprint_names:
        print("\n" + "=" * 70)
        print(f"PREDICTION | FP={fp_name} | TARGET=logS")
        print("=" * 70)

        x_data, y_data, _ = _build_xy_for_fingerprint(molecules, fp_name)
        if len(x_data) == 0:
            print(f"Aucun echantillon valide pour {fp_name}.")
            continue

        for model_name in model_names:
            stats = _evaluate_model_cv(x_data, y_data, n_splits, model_name)

            print(
                f"{model_name:>5} | "
                f"MAE={stats['mae_mean']:.3f}±{stats['mae_std']:.3f} | "
                f"R2={stats['r2_mean']:.3f}±{stats['r2_std']:.3f}"
            )

    print("\nBenchmark prédiction terminé.")
    return None
