from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.Chem.rdchem import Mol


FINGERPRINTS = {
    "morgan": lambda radius: rdFingerprintGenerator.GetMorganGenerator(radius=radius),
    "rdkit": lambda _: rdFingerprintGenerator.GetRDKitFPGenerator(),
}

SIMILARITIES = {
    "tanimoto": DataStructs.TanimotoSimilarity,
    "dice": DataStructs.DiceSimilarity,
    "cosine": DataStructs.CosineSimilarity,
}

def loadMol(mol):
    if isinstance(mol, Mol):
        return mol
    return Chem.MolFromSmiles(mol) or Chem.MolFromMolBlock(mol)


def getSimilarity(mol1, mol2, fingerprintGenerator=rdFingerprintGenerator.GetMorganGenerator(), similarity=DataStructs.TanimotoSimilarity):
    '''
    Prend un objet Mol ou le bloc de texte de la molécule.
    Retourne la similarité de Tanimoto selon la fingerprint de Morgan
    '''
    mol1, mol2 = loadMol(mol1), loadMol(mol2)
    if mol1 is None or mol2 is None:
        raise ValueError("Molécule invalide")
    
    return similarity(fingerprintGenerator.GetFingerprint(mol1), fingerprintGenerator.GetFingerprint(mol2))

#esomeprazole = Chem.MolFromSmiles('COc1ccc2nc([nH]c2c1)[S@](=O)Cc1ncc(C)c(OC)c1C')
#lansoprazole = Chem.MolFromSmiles('FC(F)(F)COc1ccnc(c1C)CS(=O)c2[nH]c3ccccc3n2')
#print(getSimilarity(esomeprazole, lansoprazole, FINGERPRINTS["morgan"](2)))