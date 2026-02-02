import pynauty
from graph import Node, Edge, MoleculeGraph
from Iso.vec_to_vc import to_vc

def to_pynauty(molGraph: MoleculeGraph, directed = False):
    # tri des ids des nodes
    nodes = list(molGraph.nodes)
    idMap = {node.id : i for i, node in enumerate(nodes)}
    n = len(nodes)

    # dict d'adjacence
    adjDict = {i : [] for i in range(n)}
    for uId, edges in molGraph.edges.items():
        if uId not in idMap: continue
        uIdx = idMap[uId]
        for edge in edges:
            vIdx = idMap[edge.v.id]
            if vIdx not in adjDict[uIdx]:
                adjDict[uIdx].append(vIdx)

    # partition de couleurs de vertex
    ColorGroups = {}
    for node in nodes:
        ColorGroups.setdefault(node.color, set()).add(idMap[node.id])
    # !!! tri des partitions de couleurs !!!
    sortedColorKeys = sorted(ColorGroups.keys())
    vertexColoring = [ColorGroups[k] for k in sortedColorKeys]

    return pynauty.Graph(n, directed, adjDict, vertexColoring)


if __name__ == "__main__":

    #test avec deux graphes id
    nodes1 = [Node(1, 'C'), Node(2, 'C'), Node(3, 'O')]
    edges1 = [Edge(nodes1[0], nodes1[1], 'double'), Edge(nodes1[1], nodes1[2], 'single')]
    g1 = MoleculeGraph(nodes1, edges1)
    tg1 = to_vc(g1, False)

    nodes2 = [Node(10, 'O'), Node(20, 'C'), Node(30, 'C')]
    edges2 = [Edge(nodes2[0], nodes2[1], 'single'), Edge(nodes2[1], nodes2[2], 'single')]
    g2 = MoleculeGraph(nodes2, edges2)

    pg1 = to_pynauty(to_vc(g1))
    pg2 = to_pynauty(to_vc(g2))

    print("Are isomorphs ? : ", pynauty.certificate(pg1) == pynauty.certificate(pg2))