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

    nbDiffAtome = 0
    nbDiffLink = 0

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
    
    def getNbDiffAtome(self):
        return self.nbDiffAtome
    
    def setNbDiffAtome(self, val):
        self.nbDiffAtome = val
    
    def getNbDiffLink(self):
        return self.nbDiffLink
    
    def setNbDiffLink(self, val):
        self.nbDiffLink = val
    
    def export_to_nauty(self):
        #à faire
        pass

    @classmethod
    def from_molfile(MoleculeGraph, path):
        nbNodes = 0
        nbEdges = 0
        find = False
        i = 1
        listNodes = []
        listEdges = []
        tempAtomes = set()
        tempLink = set()
        with open(path) as f:
            for line in f.readlines():
                splitted = line.split()
                if not find:
                    if len(splitted) < 2:
                        continue
                    if splitted[0].isdigit() and splitted[1].isdigit():
                        nbNodes, nbEdges = int(splitted[0]), int(splitted[1])
                        find = True
                else:
                    if i <= nbNodes:
                        listNodes.append(Node(i, splitted[3]))
                        tempAtomes.add(splitted[3])
                    elif i <= nbNodes + nbEdges:
                        listEdges.append(Edge(listNodes[int(splitted[0]) - 1], listNodes[int(splitted[1]) - 1], int(splitted[2])))
                        tempLink.add(splitted[2])
                    i += 1
        

        g = MoleculeGraph(listNodes, listEdges)
        g.setNbDiffAtome(len(tempAtomes))
        g.setNbDiffLink(len(tempLink))

        return g

    @classmethod
    def from_moltext(MoleculeGraph, moltext):
        nbNodes = 0
        nbEdges = 0
        find = False
        i = 1
        listNodes = []
        listEdges = []
        for line in moltext.splitlines():
            splitted = line.split()
            if not find:
                if len(splitted) < 2:
                    continue
                if splitted[0].isdigit() and splitted[1].isdigit():
                    nbNodes, nbEdges = int(splitted[0]), int(splitted[1])
                    find = True
            else:
                if i <= nbNodes:
                    listNodes.append(Node(i, splitted[3]))
                elif i <= nbNodes + nbEdges:
                    listEdges.append(Edge(listNodes[int(splitted[0]) - 1], listNodes[int(splitted[1]) - 1], splitted[2]))
                i += 1

        return MoleculeGraph(listNodes, listEdges)

if __name__ == "__main__":
    n1 = Node(1, "C")
    n2 = Node(2, "O")
    n3 = Node(3, "H")
    e1 = Edge(n1, n2, "d")
    e2 = Edge(n1, n3, "s")

    #g = MoleculeGraph([n1, n2, n3], [e1, e2])

    #print(g.get_neighbors(n1))
    #print(g.get_neighbors(n2))
    #print(g.get_neighbors(n3))
    #g = MoleculeGraph.from_molfile("C:\\Users\\ethan\\Downloads\\CHEBI_136874.mol")
    #print(g.get_neighbors(Node(1, "C")))
    import Chebi.CheBi as Chebi
    retr = Chebi.Chebi("chebi_cache.db")
    g = MoleculeGraph.from_moltext(retr.get_mol("136874"))
    print(g.get_neighbors(Node(1, "C")))
    print(g.getNbDiffAtome())
    print(g.getNbDiffLink())
