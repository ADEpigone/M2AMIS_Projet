
from distances.graph_edit_distance import GraphEditDistance

class DistanceFactory:
    """
    Factory pour créer des distances à partir de leur nom.
    """

    @staticmethod
    def create(nom_distance: str):
        if nom_distance == "test":
            return GraphEditDistance()


        raise ValueError(f"Distance inconnue : {nom_distance}")