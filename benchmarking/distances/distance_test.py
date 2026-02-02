import random
from .distance import Distance


class DistanceTest(Distance):
    """
    Distance factice pour tests : retourne une valeur pseudo-aléatoire.
    """

    def calculer_distance(self, mol1, mol2) -> float:
        # volontairement déterministe par paire pour reproductibilité
        random.seed(hash((mol1, mol2)) % 10_000)
        return random.random()
