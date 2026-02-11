import sys

import rdkit

from similarites.base.from_fingerprint_base import BaseSimilarityFromFingerprint
sys.path.append("..")
from rdkit import DataStructs
import hashlib
from collections import Counter, defaultdict
from graph import *
import math
from itertools import combinations

class CWLKernel(BaseSimilarityFromFingerprint):
    def __init__(self, iterations: int = 3, max_cycle_length: int = 8, similarity=DataStructs.TanimotoSimilarity, n_bits=2048):
        super().__init__(similarity=similarity, n_bits=n_bits)
        self.iterations = iterations
        self.max_cycle_length = max_cycle_length

    def calculate_fingerprint(self, g: MoleculeGraph) -> DataStructs:
        """
        Calcule le fingerprint à partir de l'histogramme généré par le Cellular WL.
        """
        if type(g) == DataStructs.ExplicitBitVect:
            return g
        histo = self._compute_histo(g)
        return self._counter_to_fingerprint(histo)

    def _dot(self, h1: Counter, h2: Counter) -> int:
        common_labels = set(h1.keys()) & set(h2.keys())
        score = 0
        for label in common_labels:
            score += h1[label] * h2[label]
        return score

    def _get_edge_key(self, u_id, v_id):
        """Clé canonique pour une arête (non orientée)."""
        return (min(u_id, v_id), max(u_id, v_id))

    def _build_edge_structures(self, g: MoleculeGraph):
        """
        Construit :
        - edge_labels : dict { (u_id, v_id) -> label initial de l'arête }
        - node_to_edges : dict { node_id -> set of edge_keys } (incidence)
        - edge_to_nodes : dict { edge_key -> (u_id, v_id) } (boundary)
        - edge_neighbors : dict { edge_key -> set of edge_keys } (co-boundary / upper adjacency)
        """
        edge_labels = {}
        node_to_edges = {}
        edge_to_nodes = {}

        for node in g.nodes:
            node_to_edges[node.id] = set()

        # Parcourir toutes les arêtes
        for node_id, edges in g.edges.items():
            for edge in edges:
                key = self._get_edge_key(edge.u.id, edge.v.id)
                if key not in edge_labels:
                    # Label initial d'une 1-cell : combinaison des labels des endpoints + type de liaison
                    edge_labels[key] = str(edge.color)
                    edge_to_nodes[key] = (edge.u.id, edge.v.id)
                node_to_edges.setdefault(edge.u.id, set()).add(key)
                node_to_edges.setdefault(edge.v.id, set()).add(key)

        # Co-boundary (upper adjacency) : deux arêtes sont voisines si elles partagent un nœud
        edge_neighbors = {key: set() for key in edge_labels}
        for node_id, incident_edges in node_to_edges.items():
            for e1 in incident_edges:
                for e2 in incident_edges:
                    if e1 != e2:
                        edge_neighbors[e1].add(e2)

        return edge_labels, node_to_edges, edge_to_nodes, edge_neighbors

    def _build_face_structures(self, g: MoleculeGraph, edge_labels):
        """
        Construit les structures pour les 2-cells (faces = cycles minimaux).
        - face_labels : dict { face_id -> label }
        - face_boundary : dict { face_id -> set of edge_keys } (boundary : arêtes du cycle)
        - edge_to_faces : dict { edge_key -> set of face_ids } (co-boundary des arêtes)
        - face_neighbors : dict { face_id -> set of face_ids } (adjacence : faces partageant une arête)
        """
        cycles = self._find_minimal_cycles(g)

        face_labels = {}
        face_boundary = {}
        edge_to_faces = defaultdict(set)
        face_neighbors = defaultdict(set)

        for face_id, cycle in enumerate(cycles):
            # Arêtes du cycle
            edges_in_face = set()
            edge_label_list = []
            for i in range(len(cycle)):
                ekey = self._get_edge_key(cycle[i], cycle[(i + 1) % len(cycle)])
                edges_in_face.add(ekey)
                if ekey in edge_labels:
                    edge_label_list.append(edge_labels[ekey])
                else:
                    edge_label_list.append("?")

            face_boundary[face_id] = edges_in_face

            for ekey in edges_in_face:
                edge_to_faces[ekey].add(face_id)

            # Label initial : hash des labels d'arêtes triés + taille du cycle
            edge_label_list.sort()
            init_sig = f"face_{len(cycle)}_{'_'.join(edge_label_list)}"
            face_labels[face_id] = self._hash(init_sig)

        # Adjacence entre faces : deux faces partagent au moins une arête
        for ekey, face_ids in edge_to_faces.items():
            for f1 in face_ids:
                for f2 in face_ids:
                    if f1 != f2:
                        face_neighbors[f1].add(f2)

        return face_labels, face_boundary, edge_to_faces, face_neighbors

    def _find_minimal_cycles(self, g: MoleculeGraph):
        """
        Détecte les cycles minimaux (jusqu'à max_cycle_length) dans le graphe.
        Retourne une liste de tuples de node ids (triés de manière canonique).
        """
        adj = defaultdict(set)
        for node_id, edges in g.edges.items():
            for edge in edges:
                adj[edge.u.id].add(edge.v.id)
                adj[edge.v.id].add(edge.u.id)

        found_cycles = set()

        def dfs(start, current, visited, path, depth):
            if depth > self.max_cycle_length:
                return
            for neighbor in adj[current]:
                if neighbor == start and len(path) >= 3:
                    # Cycle trouvé — canonicaliser
                    cycle = tuple(path)
                    # Rotation canonique : commencer par le plus petit id, choisir la direction minimale
                    canon = self._canonical_cycle(cycle)
                    found_cycles.add(canon)
                elif neighbor not in visited:
                    visited.add(neighbor)
                    path.append(neighbor)
                    dfs(start, neighbor, visited, path, depth + 1)
                    path.pop()
                    visited.remove(neighbor)

        for node in g.nodes:
            dfs(node.id, node.id, {node.id}, [node.id], 1)

        # Filtrer pour ne garder que les cycles minimaux (pas de cycle qui est l'union de sous-cycles)
        return self._filter_minimal_cycles(found_cycles)

    def _canonical_cycle(self, cycle):
        """Retourne la forme canonique d'un cycle (rotation min + direction min)."""
        n = len(cycle)
        rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(n)]
        reversed_cycle = tuple(reversed(cycle))
        rotations += [tuple(reversed_cycle[i:] + reversed_cycle[:i]) for i in range(n)]
        return min(rotations)

    def _filter_minimal_cycles(self, cycles):
        """Filtre les cycles pour ne garder que ceux qui ne sont pas décomposables en sous-cycles."""
        sorted_cycles = sorted(cycles, key=len)
        minimal = []
        covered_edge_sets = []

        for cycle in sorted_cycles:
            edges_of_cycle = set()
            for i in range(len(cycle)):
                e = self._get_edge_key(cycle[i], cycle[(i + 1) % len(cycle)])
                edges_of_cycle.add(e)

            # Vérifier si ce cycle est l'union exacte de cycles plus petits déjà trouvés
            is_composite = False
            if len(covered_edge_sets) >= 2:
                for r in range(2, len(covered_edge_sets) + 1):
                    if r > 3:
                        break
                    for combo in combinations(range(len(covered_edge_sets)), r):
                        union = set()
                        for idx in combo:
                            union |= covered_edge_sets[idx]
                        if union == edges_of_cycle:
                            is_composite = True
                            break
                    if is_composite:
                        break

            if not is_composite:
                minimal.append(cycle)
                covered_edge_sets.append(edges_of_cycle)

        return minimal

    def _hash(self, signature: str) -> str:
        return hashlib.sha1(signature.encode('utf-8')).hexdigest()

    def _compute_histo(self, g: MoleculeGraph) -> Counter:
        """
        Cellular WL : on raffine itérativement les labels des 0-cells (nœuds),
        1-cells (arêtes) et 2-cells (faces), puis on accumule dans un histogramme.
        """
        # Labels initiaux des 0-cells
        node_labels = {node.id: node.color for node in g.nodes}

        # Structures pour les 1-cells
        edge_labels, node_to_edges, edge_to_nodes, edge_neighbors = self._build_edge_structures(g)

        # Structures pour les 2-cells
        face_labels, face_boundary, edge_to_faces, face_neighbors = self._build_face_structures(g, edge_labels)

        # Histogramme : labels initiaux
        all_labels = []
        all_labels.extend(f"0c_{l}" for l in node_labels.values())
        all_labels.extend(f"1c_{l}" for l in edge_labels.values())
        all_labels.extend(f"2c_{l}" for l in face_labels.values())

        for _ in range(self.iterations):
            # === Mise à jour des 0-cells (nœuds) ===
            new_node_labels = {}
            for node in g.nodes:
                label_u = node_labels[node.id]

                # Adjacence (voisins via arêtes)
                adj_patterns = []
                if node.id in g.edges:
                    for edge in g.edges[node.id]:
                        neighbor_label = node_labels[edge.v.id]
                        bond_label = edge_labels[self._get_edge_key(edge.u.id, edge.v.id)]
                        adj_patterns.append(f"{neighbor_label}_{bond_label}")
                adj_patterns.sort()

                # Boundary (labels des 1-cells incidentes)
                boundary_patterns = []
                for ekey in sorted(node_to_edges.get(node.id, set())):
                    boundary_patterns.append(edge_labels[ekey])
                boundary_patterns.sort()

                signature = f"{label_u}|ADJ:{''.join(adj_patterns)}|BND:{''.join(boundary_patterns)}"
                new_node_labels[node.id] = self._hash(signature)

            # === Mise à jour des 1-cells (arêtes) ===
            new_edge_labels = {}
            for ekey, elabel in edge_labels.items():
                u_id, v_id = edge_to_nodes[ekey]

                # Boundary (labels des 2 nœuds endpoints)
                boundary = sorted([node_labels[u_id], node_labels[v_id]])

                # Co-boundary (labels des 2-cells incidentes = faces contenant cette arête)
                coboundary_patterns = []
                for fid in sorted(edge_to_faces.get(ekey, set())):
                    coboundary_patterns.append(face_labels[fid])
                coboundary_patterns.sort()

                # Adjacence (arêtes voisines partageant un nœud)
                adj_patterns = []
                for neighbor_ekey in sorted(edge_neighbors[ekey]):
                    adj_patterns.append(edge_labels[neighbor_ekey])
                adj_patterns.sort()

                signature = f"{elabel}|BND:{''.join(boundary)}|COBND:{''.join(coboundary_patterns)}|ADJ:{''.join(adj_patterns)}"
                new_edge_labels[ekey] = self._hash(signature)

            # === Mise à jour des 2-cells (faces) ===
            new_face_labels = {}
            for fid, flabel in face_labels.items():
                # Boundary (labels des arêtes du cycle)
                boundary_patterns = []
                for ekey in sorted(face_boundary[fid]):
                    boundary_patterns.append(edge_labels[ekey])
                boundary_patterns.sort()

                # Adjacence (faces partageant une arête)
                adj_patterns = []
                for neighbor_fid in sorted(face_neighbors.get(fid, set())):
                    adj_patterns.append(face_labels[neighbor_fid])
                adj_patterns.sort()

                signature = f"{flabel}|BND:{''.join(boundary_patterns)}|ADJ:{''.join(adj_patterns)}"
                new_face_labels[fid] = self._hash(signature)

            # Mettre à jour tous les labels
            node_labels = new_node_labels
            edge_labels = new_edge_labels
            face_labels = new_face_labels

            # Accumuler dans l'histogramme
            all_labels.extend(f"0c_{l}" for l in node_labels.values())
            all_labels.extend(f"1c_{l}" for l in edge_labels.values())
            all_labels.extend(f"2c_{l}" for l in face_labels.values())

        return Counter(all_labels)


if __name__ == "__main__":

    # Exemple avec un cycle (triangle) pour tester les 2-cells
    n1 = Node(1, "C")
    n2 = Node(2, "C")
    n3 = Node(3, "O")
    g1 = MoleculeGraph([n1, n2, n3], [Edge(n1, n2, "1"), Edge(n2, n3, "1"), Edge(n1, n3, "1")])

    n4 = Node(1, "C")
    n5 = Node(2, "C")
    n6 = Node(3, "O")
    g2 = MoleculeGraph([n4, n5, n6], [Edge(n4, n5, "2"), Edge(n5, n6, "1"), Edge(n4, n6, "1")])

    n7 = Node(10, "C")
    n8 = Node(11, "C")
    n9 = Node(12, "O")
    g3 = MoleculeGraph([n7, n8, n9], [Edge(n7, n8, "1"), Edge(n8, n9, "1"), Edge(n7, n9, "1")])

    cwl = CWLKernel(iterations=2)
    print(cwl._compute_histo(g1))
    sim_1_2 = cwl.calculate_similarity(g1, g2)
    sim_1_3 = cwl.calculate_similarity(g1, g3)

    print(f"Similarité CWL G1 G2 : {sim_1_2}")
    print(f"Similarité CWL G1 G3 : {sim_1_3}")