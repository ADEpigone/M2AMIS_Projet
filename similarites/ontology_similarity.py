import math
from typing import Optional, Set
from Chebi.ontology.ontology_tree import OntologyNode, OntologyTree
from graph import MoleculeGraph
from similarites.base.graph_similarity import BaseGraphSimilarity

class OntologySimilarity(BaseGraphSimilarity):
    
    def __init__(self, ontology, decay_factor=0.1):
        self.ontology = ontology
        self.decay_factor = decay_factor
        
    def calculate_similarity(self, g1, g2) -> float:
        id_a = g1.chebi_id
        id_b = g2.chebi_id
        
        if id_a == id_b:
            return 1.0
            
        node_a = self.ontology.get_node(id_a)
        node_b = self.ontology.get_node(id_b)
        if not node_a or not node_b:
            return 0.0

        lca = self.ontology.get_lca(id_a, id_b)
        if not lca:
            return 0.0
            
        d_a = node_a.get_depth()
        d_b = node_b.get_depth()
        d_lca = lca.get_depth()
        
        distance = (d_a + d_b) - (2 * d_lca)
        
        sim = math.exp(-self.decay_factor * distance)
        
        return sim