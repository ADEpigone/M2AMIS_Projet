from .distance_test import DistanceTest


class DistanceFactory:
    """
    Factory pour créer des distances à partir de leur nom.
    """

    @staticmethod
    def create(nom_distance: str):
        if nom_distance == "test":
            return DistanceTest()
        

        raise ValueError(f"Distance inconnue : {nom_distance}")
