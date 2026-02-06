import sys

from similarites.base.graph_similarity import BaseGraphSimilarity
sys.path.append("..")

import hashlib
from collections import Counter
from graph import *  # Assure-toi que l'import correspond à ton fichier
import math
class WLKernel(BaseGraphSimilarity):
    def __init__(self, iterations: int = 3):
        self.iterations = iterations

    def calculate_similarity(self, g1: MoleculeGraph, g2: MoleculeGraph) -> int:
        """
        Calcule la similarité entre deux graphes moléculaires g1 et g2
        à l'aide 1-WL t'as capté
        """
        hist_g1 = self._compute_histo(g1)
        hist_g2 = self._compute_histo(g2)

        return self._dot(hist_g1, hist_g2)/math.sqrt(self._dot(hist_g1, hist_g1) * self._dot(hist_g2, hist_g2))

    def _dot(self, h1: Counter, h2: Counter) -> int:
        """
        Produit scalaire entre 2 histogrammes
        On compare que les labels en commun
        Là ceux qui sont pas en commun servent à rien c'est un peu l'intersection pour la sim finale
        """

        common_labels = set(h1.keys()) & set(h2.keys())
        
        score = 0
        for label in common_labels:
            score += h1[label] * h2[label]
            
        return score

    def _compute_histo(self, g: MoleculeGraph) -> list[str]:
        """
        Génère la liste de TOUS les labels rencontrés pdt iterations étapes.
        """
        current_labels = {node.id: node.color for node in g.nodes}
        
        all_labels = list(current_labels.values())

        for _ in range(self.iterations):
            new_labels = {}
            
            for node in g.nodes:
                label_u = current_labels[node.id]
                
                neighbors_pattern = []
                
                #on ajoute les labels des voisins avec le type de liaison (label sommet + _ + label liaison, on mais jamais à jour les liaisons)
                if node.id in g.edges:
                    for edge in g.edges[node.id]:
                        neighbor_label = current_labels[edge.v.id]
                        bond_type = str(edge.color)
                        neighbors_pattern.append(f"{neighbor_label}_{bond_type}")
                
                neighbors_pattern.sort()
                
                signature = label_u + ":" + "".join(neighbors_pattern)
                
                #on hash, hash() pas ultra fiable donc sha1 ou sha256 ou whatever tbh
                hashed_label = hashlib.sha1(signature.encode('utf-8')).hexdigest()
                new_labels[node.id] = hashed_label
            
            current_labels = new_labels
            all_labels.extend(current_labels.values())
            
        return Counter(all_labels)

if __name__ == "__main__":
    
    n1 = Node(1, "C")
    n2 = Node(2, "O")
    g1 = MoleculeGraph([n1, n2], [Edge(n1, n2, "1")]) 

    n3 = Node(1, "C")
    n4 = Node(2, "O")
    g2 = MoleculeGraph([n3, n4], [Edge(n3, n4, "2")])

    n5 = Node(10, "C")
    n6 = Node(11, "O")
    g3 = MoleculeGraph([n5, n6], [Edge(n5, n6, "1")])

    wl = WLKernel(iterations=2)

    sim_1_2 = wl.calculate_similarity(g1, g2)
    sim_1_3 = wl.calculate_similarity(g1, g3)

    print(f"Similarité G1 G2 : {sim_1_2}") 
    print(f"Similarité G1 G3: {sim_1_3}")