from typing import Dict, List
from .distances.distance import Distance


def benchmark_knn(
    distance: Distance,
    dataset: Dict[str, List[str]],
    k: int = 3
):
    """
    Benchmark kNN avec vérité terrain figée.

    dataset = {
        "mol1": ["voisin1", "voisin2", ...],
        ...
    }
    """
    resultats = {}

    for mol, voisins_attendus in dataset.items():
        distances = []
        for autre in dataset.keys():
            if autre == mol:
                continue
            d = distance.calculer_distance(mol, autre)
            distances.append((autre, d))

        distances.sort(key=lambda x: x[1])
        voisins_predits = [m for m, _ in distances[:k]]

        intersection = set(voisins_predits) & set(voisins_attendus)

        precision_k = len(intersection) / k
        
        if len(voisins_attendus) > 0:
            rappel_k = len(intersection) / len(voisins_attendus)
        else:
            rappel_k = None  # ou float("nan")


        resultats[mol] = {
            "precision@k": precision_k,
            "rappel@k": rappel_k,
            "voisins_predits": voisins_predits,
        }

    return resultats


def benchmark_monotonie_separation(
    distance: Distance,
    dataset: List[str]
):
    """
    Dataset attendu : [seed, e1, e2, e3, random]
    """
    seed, e1, e2, e3, random_mol = dataset

    d1 = distance.calculer_distance(seed, e1)
    d2 = distance.calculer_distance(seed, e2)
    d3 = distance.calculer_distance(seed, e3)
    dr = distance.calculer_distance(seed, random_mol)

    monotonie = d1 < d2 < d3
    separation = max(d1, d2, d3) < dr

    return {
        "distances": {
            "d(seed,e1)": d1,
            "d(seed,e2)": d2,
            "d(seed,e3)": d3,
            "d(seed,random)": dr,
        },
        "monotonie_ok": monotonie,
        "separation_ok": separation,
    }
