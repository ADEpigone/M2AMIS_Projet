import sys
sys.path.append("..")
from distances.base.graph_distance import BaseGraphDistance
from graph import MoleculeGraph

class GraphEditDistance(BaseGraphDistance):
    def calculate_distance(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        if len(g1.nodes) > len(g2.nodes):
            g1, g2 = g2, g1

        nodes_g1 = list(g1.nodes)
        nodes_g2 = list(g2.nodes)
        
        return self._rec(0, {}, nodes_g1, nodes_g2, g1, g2)

    def _rec(self, idx_g1, mapping, nodes_g1, nodes_g2, g1, g2):
        """
        idx_g1  : L'index de l'atome de G1 que l'on traite
        mapping : Un dictionnaire {index_g1: index_g2} des choix faits jusqu'ici
                  si on a i : None -> i a été suppr
        """
        
        if idx_g1 == len(nodes_g1):
            return self._calculate_remaining_cost(mapping, nodes_g2, g2)
        
        u = nodes_g1[idx_g1]
        best_cost = float('inf')

        used_g2_indices = set(v_idx for v_idx in mapping.values() if v_idx is not None)

        # on regarde tous les u -> v possibles (on map)
        for idx_g2, v in enumerate(nodes_g2):
            if idx_g2 not in used_g2_indices:
                
                #si on doit changer de couleur/type atome
                node_cost = 0 if u.color == v.color else 1
                
                edge_cost = self._calculate_edge_cost(idx_g1, idx_g2, mapping, g1, g2, nodes_g1, nodes_g2)
                
                current_step_cost = node_cost + edge_cost
                
                #on register le mapping et on rec dessus
                mapping[idx_g1] = idx_g2
                total_cost = current_step_cost + self._rec(idx_g1 + 1, mapping, nodes_g1, nodes_g2, g1, g2)
                
                #on annule, c du backtrack quoi
                del mapping[idx_g1]
                
                best_cost = min(best_cost, total_cost)

        #ici on regarde juste ce qu'il se passe quand on suppr
        #on backtrack et gg
        #j'ai la flm de continuer els commentaires
        del_edge_cost = self._calculate_edge_del_cost(idx_g1, mapping, g1, nodes_g1)
        
        mapping[idx_g1] = None
        total_cost_del = 1 + del_edge_cost + self._rec(idx_g1 + 1, mapping, nodes_g1, nodes_g2, g1, g2)
        del mapping[idx_g1]

        best_cost = min(best_cost, total_cost_del)

        return best_cost

    def _calculate_edge_cost(self, u_idx, v_idx, mapping, g1, g2, nodes_g1, nodes_g2):
        """
        Etant donné u -> v, calcule le coûts des arêtes déjà mappées
        <-> cohérence topologique, on regarde pour notre nouveau mapping si les liens avec u et v sont cohérents (si non on punit)
        """
        cost = 0
        u = nodes_g1[u_idx]
        v = nodes_g2[v_idx]

        for prev_u_idx, prev_v_idx in mapping.items():
            prev_u = nodes_g1[prev_u_idx]
            

            if prev_v_idx is None:
                if g1.has_edge(u, prev_u):
                    cost += 1 
                continue

            prev_v = nodes_g2[prev_v_idx]
            
            has_edge_g1 = g1.has_edge(u, prev_u)
            has_edge_g2 = g2.has_edge(v, prev_v)

            #si on a pas la même structure
            #ou si jamais on a la même structure (mais pas le même type et qu'on doit map dcp)
            if has_edge_g1 != has_edge_g2:
                cost += 1
            elif has_edge_g1 and has_edge_g2:
                bond_g1 = g1.get_edge_data(u, prev_u)
                bond_g2 = g2.get_edge_data(v, prev_v)
                if bond_g1 != bond_g2:
                    cost += 1

        return cost

    def _calculate_edge_del_cost(self, u_idx, mapping, g1, nodes_g1):
        """
        Etant donné un sommet qu'on veut suppr on va regarder cmb d'arêtes dans le mapping on doit suppr
        """
        cost = 0
        u = nodes_g1[u_idx]
        for prev_u_idx, _ in mapping.items():
            prev_u = nodes_g1[prev_u_idx]
            if g1.has_edge(u, prev_u):
                cost += 1
        return cost

    def _calculate_remaining_cost(self, mapping, nodes_g2, g2):
        """
        On calcule juste le surplus, les choses qu'on a pas pu map/prendre en compte
        """
        cost = 0
        used_g2 = set(v for v in mapping.values() if v is not None)
        
        remaining_indices = [i for i in range(len(nodes_g2)) if i not in used_g2]
        cost += len(remaining_indices)
        
        for i in remaining_indices:
            u_rest = nodes_g2[i]
            
            for j in remaining_indices:
                if j > i:
                    v_rest = nodes_g2[j]
                    if g2.has_edge(u_rest, v_rest):
                        cost += 1
            
            for k in used_g2:
                v_used = nodes_g2[k]
                if g2.has_edge(u_rest, v_used):
                    cost += 1
                    
        return cost

if __name__ == "__main__":
    import sys 
    sys.path.append("../")
    from Chebi.CheBi import Chebi
    from graph import MoleculeGraph
    
    
    retr = Chebi("chebi_cache.db")
    g1 = MoleculeGraph.from_moltext(retr.get_mol("15377"))
    g2 = MoleculeGraph.from_moltext(retr.get_mol("15378"))

    ged = GraphEditDistance()
    print("Distance:", ged.calculate_distance(g1, g2))