from graph import Node, Edge, MoleculeGraph

def transformation(graph : MoleculeGraph, debug = False):

    # récup des couleurs
    colorsTemp = set()
    k, values = zip(*graph.edges.items())
    print(list(values))

    for edges in list(values):
        for edge in edges:
            colorsTemp.add(edge.color)

    # tri des couleurs et mapping vers int
    colors = sorted(list(colorsTemp))
    colorsToInt = {color : i + 1 for i, color in enumerate(colors)}
    print("colors : " + str(colors))
    print("colortoint : " + str(colorsToInt))

    nbOfColors = len(colors)
    if nbOfColors == 0:
        nbOfLayers = 1
    else:
        nbOfLayers = nbOfColors.bit_length()
    print("Nb layers : " + str(nbOfLayers))

    N = len(graph.nodes)
    newNodes = []

    # helper pour ids des nouveaux nodes
    def get_new_id(ogId, layerId):
        return (layerId * N) + ogId
    
    # nouveaux nodes
    for layerId in range(nbOfLayers):
        for node in graph.nodes:
            newId = get_new_id(node.id, layerId)
            newColor = f"{node.color}_layer{layerId}"
            newNodes.append(Node(newId, newColor))

    # nouveaux edges
    nodeMap = {n.id : n for n in newNodes}
    newEdges = []

    # edges verticaux
    for layerId in range(nbOfLayers - 1):
        for node in graph.nodes:
            u = nodeMap[get_new_id(node.id, layerId)]
            v = nodeMap[get_new_id(node.id, layerId + 1)]
            newEdges.append(Edge(u, v, "verticalLink"))

    
    # edges horizontaux
    for edges in list(values):
        print(edges)
        for edge in edges:
            colorVal = colorsToInt[edge.color]

            for bit in range(nbOfLayers):
                if (colorVal >> bit) & 1:
                    uNew = nodeMap[get_new_id(edge.u.id, bit)]
                    vNew = nodeMap[get_new_id(edge.v.id, bit)]
                    newEdges.append(Edge(uNew, vNew, "horizontalLink"))
    
    if debug == True:
        print("NODES : ")
        for node in newNodes:
            print("N: ", str(node))
        print("#############################################")
        print("EDGES : ")
        for edge in newEdges:
            print("E: ", str(edge))
    
    return MoleculeGraph(newNodes, newEdges)

def test_node_and_colors(graph1, graph2, nboflayer):
    print("Len nodes : " , len(graph2.nodes) == len(graph1.nodes) * nboflayer)

    nbOfCol1 = set()
    for node in graph1.nodes:
        nbOfCol1.add(node.color)
    nbOfCol2 = set()
    for node in graph2.nodes:
        nbOfCol2.add(node.color)
    print("Nb of colors : " , len(nbOfCol2) == (len(nbOfCol1) * nboflayer))
    

if __name__ == "__main__":
    g = MoleculeGraph.molFromFile("C:\\Users\\cyria\\Downloads\\CHEBI_27518.mol")
    
    newG = transformation(g, True)
    test_node_and_colors(g, newG, 1)