from similarites.base.factory import DistanceFactory
from .methodes import (
    benchmark_knn,
    benchmark_monotonie_separation,
)


def main_benchmark():
    # Création de la distance via la factory
    distance = DistanceFactory.create("test")

    # Dataset figé pour kNN
    dataset_knn = {
        "mol_A": ["mol_B", "mol_C"],
        "mol_B": ["mol_A"],
        "mol_C": ["mol_A"],
        "mol_D": [],
    }

    print("=== Benchmark kNN ===")
    resultats_knn = benchmark_knn(distance, dataset_knn, k=2)
    for mol, res in resultats_knn.items():
        print(mol, res)

    # Dataset figé pour monotonie / séparation
    dataset_monotonie = ["seed", "e1", "e2", "e3", "random"]

    print("\n=== Benchmark Monotonie / Séparation ===")
    res_mono = benchmark_monotonie_separation(distance, dataset_monotonie)
    print(res_mono)


if __name__ == "__main__":
    main_benchmark()
