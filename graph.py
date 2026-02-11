from typing import NamedTuple
from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

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

    def __init__(self, nodes : list[Node] = None, edges : list[Edge] = None, directed=False, mol_file = None, chebi_id = None):
        self.nodes = set()
        self.mol = mol_file
        self.chebi_id = chebi_id
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

    def has_edge(self, u : Node, v : Node) -> bool:
        if u.id not in self.edges:
            return False
        for edge in self.edges[u.id]:
            if edge.v == v:
                return True
        return False

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

    @staticmethod
    @lru_cache(maxsize=1024)
    def _rdkit_props_from_molblock(molblock: str) -> dict:
        # Cache calculs RDKit a partir du molblock.
        if not molblock:
            return {}
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        if mol is None:
            return {}
        return {
            "logP": float(Crippen.MolLogP(mol)),
            "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
            "mol_wt": float(Descriptors.MolWt(mol)),
        }

    @staticmethod
    @lru_cache(maxsize=2048)
    def _molblock_from_smiles(smiles: str) -> str:
        if not smiles:
            return ""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        return Chem.MolToMolBlock(mol)

    def get_properties(self) -> dict:
        if not self.mol:
            return {}
        return self._rdkit_props_from_molblock(self.mol)

    @staticmethod
    def _get_atom_label_rich(atom) -> str:
        """
        Symbole + H + Aromaticité.
        
        Ex: "C_h1_a1" (au lieu de "C_d3_h1_a1_c0")
        """
        symbol = atom.GetSymbol()
        
        # Le degré est implicite dans le graphe, on le retire du label dur.
        # degree = atom.GetTotalDegree() 
        
        # Les H sont cruciaux pour les donneurs/accepteurs (OH vs =O)
        hs = atom.GetTotalNumHs() 
        
        # L'aromaticité est cruciale pour la rigidité/solubilité
        is_arom = "1" if atom.GetIsAromatic() else "0"
        
        # La charge est souvent gérée par le contexte (N+ a 4 voisins)
        # charge = atom.GetFormalCharge()
        
        return f"{symbol}_h{hs}_a{is_arom}"

    @classmethod
    def from_moltext(cls, moltext: str, chebi_id=None):
        """
        Parse le bloc MOL avec RDKit pour extraire une topologie riche.
        """
        if not moltext:
            return cls([], [], mol_file=moltext, chebi_id=chebi_id)

        # 1. Parsing RDKit
        mol = Chem.MolFromMolBlock(moltext, sanitize=True, removeHs=False)
        
        # Si le molblock est pourri, fallback ou erreur (ici on renvoie vide)
        if mol is None:
            # Optionnel: Essayer sans sanitize si ça plante souvent
             mol = Chem.MolFromMolBlock(moltext, sanitize=False, removeHs=False)
        
        if mol is None:
            return cls([], [], mol_file=moltext, chebi_id=chebi_id)

        # 2. Création des Nodes Enrichis
        rdkit_nodes = {} # Map idx RDKit -> Node object
        node_list = []
        
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            # C'EST LA QUE LA MAGIE OPERE : LABEL RICHE
            label = cls._get_atom_label_rich(atom) 
            
            # Note: RDKit idx commence à 0, tes ids précédents commençaient à 1
            # On garde 0-indexed c'est plus standard, ou atom.GetIdx()+1 si tu veux garder la compatibilité legacy
            n = Node(idx, label)
            rdkit_nodes[idx] = n
            node_list.append(n)

        # 3. Création des Edges
        edge_list = []
        for bond in mol.GetBonds():
            idx_u = bond.GetBeginAtomIdx()
            idx_v = bond.GetEndAtomIdx()
            
            u = rdkit_nodes[idx_u]
            v = rdkit_nodes[idx_v]
            
            # Enrichissement possible du lien aussi (Single/Double/Aromatic)
            btype = str(bond.GetBondType())
            
            edge_list.append(Edge(u, v, btype))

        return cls(node_list, edge_list, mol_file=moltext, chebi_id=chebi_id)

    @classmethod
    def from_molfile(cls, path: str, chebi_id=None):
        with open(path, "r") as f:
            moltext = f.read()
        return cls.from_moltext(moltext, chebi_id=chebi_id)

    @classmethod
    def from_smiles(cls, smiles: str, chebi_id=None):
        moltext = cls._molblock_from_smiles(smiles)
        return cls.from_moltext(moltext, chebi_id=chebi_id)

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
