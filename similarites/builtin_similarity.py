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
    def __init__(self, fingerprint = "morgan", similarity ="tanimoto"):
        if type(fingerprint) == str and fingerprint in FINGERPRINTS:
            self.fingerprint = FINGERPRINTS[fingerprint](2)
        else:
            raise ValueError(f"Fingerprint {fingerprint} non supportée. Choisissez parmi : {list(FINGERPRINTS.keys())}")
        super().__init__(similarity)

    def calculate_fingerprint(self, g1: MoleculeGraph) -> DataStructs:
        if type(g1) == DataStructs.ExplicitBitVect:
            return g1
        if g1.mol is None:
            raise ValueError("Molécule sans mol stocké, que faire ?")
        mol = g1.mol
        m = Chem.MolFromMolBlock(mol, sanitize=True)

        return self.fingerprint.GetFingerprint(m)

#BAH DU COUP CECI MARCHE PAS : IL FAUT PASSER PAR DB QUAND EWAN AURA MERGE
# PUNAISE MERGE SRX
# mais l'use est la même j'ai juste encapsulé

#esomeprazole = Chem.MolFromSmiles('COc1ccc2nc([nH]c2c1)[S@](=O)Cc1ncc(C)c(OC)c1C')
#lansoprazole = Chem.MolFromSmiles('FC(F)(F)COc1ccnc(c1C)CS(=O)c2[nH]c3ccccc3n2')
#print(getSimilarity(esomeprazole, lansoprazole, FINGERPRINTS["morgan"](2)))