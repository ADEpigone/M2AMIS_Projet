from similarites.base.from_fingerprint_base import BaseSimilarityFromFingerprint

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.Chem.rdchem import Mol
from graph import MoleculeGraph

FINGERPRINTS = {
    "morgan": lambda radius: rdFingerprintGenerator.GetMorganGenerator(radius=radius),
    "rdkit": lambda _: rdFingerprintGenerator.GetRDKitFPGenerator(),
}


SIMILARITIES = {
    "tanimoto": DataStructs.TanimotoSimilarity,
    "dice": DataStructs.DiceSimilarity,
    "cosine": DataStructs.CosineSimilarity,
}


class BuiltinSimilarity(BaseSimilarityFromFingerprint):
    '''
    Prend un objet Mol ou le bloc de texte de la molécule.
    Retourne la similarité de Tanimoto selon la fingerprint de Morgan
    '''
    def __init__(self, fingerprint_function = "morgan", similarity_function="tanimoto"):
        if type(fingerprint_function) == str and fingerprint_function in FINGERPRINTS:
            self.fingerprint_function = FINGERPRINTS[fingerprint_function]
        if type(similarity_function) == str and similarity_function in SIMILARITIES:
            similarity_function = SIMILARITIES[similarity_function]
        super().__init__(similarity_function)

    def calculate_fingerprint(self, g1: MoleculeGraph, g2: MoleculeGraph) -> float:
        if g1.mol is None or g2.mol is None:
            raise ValueError("Molécule sans mol stocké, que faire ?")
        mol1, mol2 = g1.mol, g2.mol

        return self.fingerprint_function.GetFingerprint(mol1), self.fingerprint_function.GetFingerprint(mol2)
    

#BAH DU COUP CECI MARCHE PAS : IL FAUT PASSER PAR DB QUAND EWAN AURA MERGE
# PUNAISE MERGE SRX
# mais l'use est la même j'ai juste encapsulé

#esomeprazole = Chem.MolFromSmiles('COc1ccc2nc([nH]c2c1)[S@](=O)Cc1ncc(C)c(OC)c1C')
#lansoprazole = Chem.MolFromSmiles('FC(F)(F)COc1ccnc(c1C)CS(=O)c2[nH]c3ccccc3n2')
#print(getSimilarity(esomeprazole, lansoprazole, FINGERPRINTS["morgan"](2)))