from typing import NamedTuple

# https://benhoyt.com/writings/python-pattern-matching/
# énorme masterclass j'avais envie de tenter


class Node(NamedTuple):
    id: int
    color: str


class Edge(NamedTuple):
    u: Node
    v: Node
    color : str

class MoleculeGraph:
    def __init__(self, nodes : list[Node] = None, edges : list[Edge] = None, directed=False):
        self.nodes = set()
        self.edges = {}
        #on pourrait utiliser implicitement les nodes à partir de ceux des edges, mais je pense que pour l'instant c'est mieux
        for elem in nodes:
            self.add_node(elem)
        for elem in edges:
            self.add_edge(elem, directed=directed)
        

    def add_node(self, node : Node):
        self.nodes.add(node)

    def add_edge(self, edge : Edge, directed=False):
        if edge.u.id not in self.edges:
            self.edges[edge.u.id] = []
        if not directed and edge.v.id not in self.edges:
            self.edges[edge.v.id] = []
        self.edges[edge.u.id].append(edge)
        if not directed:
            self.edges[edge.v.id].append(Edge(edge.v, edge.u, edge.color))

    def get_neighbors(self, node : Node) -> list[Node]:
        if node.id not in self.edges:
            return []
        return [edge.v for edge in self.edges[node.id]]
    
    def export_to_nauty(self):
        #à faire
        pass

if __name__ == "__main__":
    n1 = Node(1, "C")
    n2 = Node(2, "O")
    n3 = Node(3, "H")
    e1 = Edge(n1, n2, "d")
    e2 = Edge(n1, n3, "s")

    g = MoleculeGraph([n1, n2, n3], [e1, e2])

    print(g.get_neighbors(n1))
    print(g.get_neighbors(n2))
    print(g.get_neighbors(n3))