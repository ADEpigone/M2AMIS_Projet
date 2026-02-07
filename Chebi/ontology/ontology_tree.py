from __future__ import annotations
from typing import Optional


class OntologyNode:
    """Nœud de l'arbre d'ontologie ChEBI."""

    def __init__(self, chebi_id: str, name: str = "", definition: str = ""):
        self.chebi_id = chebi_id
        self.name = name
        self.definition = definition
        self.synonyms: list[str] = []
        self.properties: dict[str, list[str]] = {}  # relationship type -> list of target ids
        self.parents: list[OntologyNode] = []  # is_a parents
        self.children: list[OntologyNode] = []

    def add_parent(self, parent: OntologyNode):
        if parent not in self.parents:
            self.parents.append(parent)
        if self not in parent.children:
            parent.children.append(self)

    def add_relationship(self, rel_type: str, target_id: str):
        if rel_type not in self.properties:
            self.properties[rel_type] = []
        self.properties[rel_type].append(target_id)

    def get_ancestors(self) -> set[str]:
        """Retourne tous les ancêtres (transitivement) via is_a."""
        ancestors = set()
        stack = list(self.parents)
        while stack:
            node = stack.pop()
            if node.chebi_id not in ancestors:
                ancestors.add(node.chebi_id)
                stack.extend(node.parents)
        return ancestors

    def get_descendants(self) -> set[str]:
        """Retourne tous les descendants (transitivement)."""
        descendants = set()
        stack = list(self.children)
        while stack:
            node = stack.pop()
            if node.chebi_id not in descendants:
                descendants.add(node.chebi_id)
                stack.extend(node.children)
        return descendants

    def get_depth(self) -> int:
        """Profondeur minimale depuis la racine."""
        if not self.parents:
            return 0
        return 1 + min(p.get_depth() for p in self.parents)

    def __repr__(self):
        return f"OntologyNode({self.chebi_id}, {self.name!r})"


class OntologyTree:
    """Arbre d'ontologie ChEBI complet."""

    def __init__(self):
        self.nodes: dict[str, OntologyNode] = {}
        self.roots: list[OntologyNode] = []

    def get_or_create_node(self, chebi_id: str) -> OntologyNode:
        if chebi_id not in self.nodes:
            self.nodes[chebi_id] = OntologyNode(chebi_id)
        return self.nodes[chebi_id]

    def get_node(self, chebi_id: str) -> Optional[OntologyNode]:
        return self.nodes.get(chebi_id)

    def build_roots(self):
        """Identifie les racines (nœuds sans parents)."""
        self.roots = [n for n in self.nodes.values() if not n.parents]

    def get_lca(self, id_a: str, id_b: str) -> Optional[OntologyNode]:
        """Lowest Common Ancestor entre deux nœuds."""
        node_a = self.get_node(id_a)
        node_b = self.get_node(id_b)
        if not node_a or not node_b:
            return None

        ancestors_a = node_a.get_ancestors()
        ancestors_a.add(id_a)
        ancestors_b = node_b.get_ancestors()
        ancestors_b.add(id_b)

        common = ancestors_a & ancestors_b
        if not common:
            return None

        # Le LCA est l'ancêtre commun le plus profond
        best = None
        best_depth = -1
        for cid in common:
            node = self.get_node(cid)
            if node:
                d = node.get_depth()
                if d > best_depth:
                    best_depth = d
                    best = node
        return best

    def search_by_name(self, query: str) -> list[OntologyNode]:
        """Recherche par nom (case-insensitive, substring)."""
        query_lower = query.lower()
        return [n for n in self.nodes.values() if query_lower in n.name.lower()]

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return f"OntologyTree({len(self.nodes)} nodes, {len(self.roots)} roots)"