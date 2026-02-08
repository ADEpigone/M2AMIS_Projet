import numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import squareform
from tqdm import tqdm
import json
import random
from Chebi.CheBi2 import CheBi2
from utils import load_ontology
from similarites.cwl_kernel import CWLKernel
from similarites.ontology_similarity import OntologySimilarity
from clustering.clusters import *
from graph import MoleculeGraph
from similarites.builtin_similarity import BuiltinSimilarity
def compute_correlation_syntactic_semantic(molecules, syn_kernel, ontology_sim):
    """
    Calcule la corrélation entre similarités syntaxiques et sémantiques.
    
    Returns:
        dict avec corrélations Pearson, Spearman, et scatter data
    """
    print("\n" + "="*80)
    print("SYNTAXIQUE vs SÉMANTIQUE")
    print("="*80)
    
    n = len(molecules)
    print(f"\nCalcul des similarités pour {n} molécules ({n*(n-1)//2} paires)...")
    
    print("Calcul des fingerprints ...")
    fingerprints = [syn_kernel.calculate_fingerprint(m["graph"]) for m in tqdm(molecules, desc="Fingerprints")]
    
    syntactic_sims = []
    semantic_sims = []
    
    print("Calcul des similarités...")
    for i in tqdm(range(n), desc="Avancement"):
        for j in range(i + 1, n):
            sim_syn = syn_kernel.calculate_similarity(fingerprints[i], fingerprints[j])
            syntactic_sims.append(sim_syn)
            
            sim_sem = ontology_sim.calculate_similarity(molecules[i]["graph"], molecules[j]["graph"])
            semantic_sims.append(sim_sem)
    
    syntactic_sims = np.array(syntactic_sims)
    semantic_sims = np.array(semantic_sims)
    
    pearson_r, pearson_p = pearsonr(syntactic_sims, semantic_sims)
    spearman_r, spearman_p = spearmanr(syntactic_sims, semantic_sims)
    
    results = {
        "n_molecules": n,
        "n_pairs": len(syntactic_sims),
        "pearson_r": float(pearson_r),
        "spearman_r": float(spearman_r),
        "syntactic_mean": float(syntactic_sims.mean()),
        "syntactic_std": float(syntactic_sims.std()),
        "semantic_mean": float(semantic_sims.mean()),
        "semantic_std": float(semantic_sims.std())
    }
    
    print("\nRÉSULTATS:")
    print(f"  Corrélation de Pearson:   r = {pearson_r:.4f} (p = {pearson_p:.2e})")
    print(f"  Corrélation de Spearman:  ρ = {spearman_r:.4f} (p = {spearman_p:.2e})")
    
    print(f"\nSTATISTIQUES:")
    print(f"  Similarité syntaxique:   {syntactic_sims.mean():.3f} ± {syntactic_sims.std():.3f}")
    print(f"  Similarité sémantique:   {semantic_sims.mean():.3f} ± {semantic_sims.std():.3f}")
    
    print(f"\nCORRÉLATION PAR NIVEAUX DE SIMILARITÉ SYNTAXIQUE:")
    bins = [(0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
    for low, high in bins:
        mask = (syntactic_sims >= low) & (syntactic_sims < high)
        if mask.sum() > 0:
            avg_sem = semantic_sims[mask].mean()
            count = mask.sum()
            print(f"  Syntaxique [{low:.1f}-{high:.1f}]: {count:6d} paires / Sémantique moyen = {avg_sem:.3f}")
    
    print("="*80 + "\n")
    
    return results

if __name__ == "__main__":

    db = CheBi2(DB_PATH)
    ontology = load_ontology()
    
    #MERCI GEMINI !!!!
    TARGET_CLASSES = {
        # --- GROUPE 1: LES CLASSIQUES (Biochimie) ---
        "Amino acid": "CHEBI:33709",       # Squelette N-C-C
        "Fatty acid": "CHEBI:35366",       # Chaîne aliphatique
        "Monosaccharide": "CHEBI:35381",   # Polyol cyclique
        "Steroid": "CHEBI:35341",          # 4 cycles fusionnés
        "Nucleobase": "CHEBI:18059",       # Bases A, T, C, G, U
        "Phospholipid": "CHEBI:16247",     # Tête polaire + queues grasses
        "Vitamin D": "CHEBI:27300",        # Secostéroides (cycles ouverts)
        "Prostaglandin": "CHEBI:26333",    # Cycle à 5 + chaînes
        "Sphingolipid": "CHEBI:26739",     # Squelette sphingosine
        "Carotenoid": "CHEBI:23044",       # Longue chaîne conjuguée

        # --- GROUPE 2: ANTIBIOTIQUES & PHARMA (Scaffolds complexes) ---
        "Penicillin": "CHEBI:17234",       # Beta-lactame + cycle à 5
        "Cephalosporin": "CHEBI:23066",    # Beta-lactame + cycle à 6
        "Tetracycline": "CHEBI:26895",     # 4 cycles linéaires
        "Macrolide": "CHEBI:25106",        # Lactone macrocyclique (géant)
        "Sulfonamide": "CHEBI:35358",      # Groupe S(=O)2-N
        "Benzodiazepine": "CHEBI:22720",   # Bicyclique diazépine (psychoactif)
        "Barbiturate": "CHEBI:22693",      # Cycle pyrimidine-trione
        "Quinolone": "CHEBI:26523",        # Bicyclique 4-quinolone
        "Morphinan": "CHEBI:36362",        # Squelette opioïde (ponté)
        "Phenothiazine": "CHEBI:37943",    # Tricyclique S-N (antipsychotique)

        # --- GROUPE 3: NATUREL & PLANTES (Cycles aromatiques) ---
        "Flavonoid": "CHEBI:47916",        # C6-C3-C6
        "Coumarin": "CHEBI:23377",         # Benzopyrone
        "Isoflavonoid": "CHEBI:50753",     # Variante flavonoides
        "Anthocyanin": "CHEBI:38697",      # Pigments colorés
        "Catechin": "CHEBI:23053",         # Flavan-3-ols
        "Lignan": "CHEBI:25030",           # Dimères de phénylpropanoïdes
        "Stilbene": "CHEBI:26780",         # Diarylethène (ex: Resveratrol)
        "Alkaloid": "CHEBI:22315",         # Contient de l'Azote (très varié)
        "Tannin": "CHEBI:26848",           # Polyphénols complexes

        # --- GROUPE 4: TERPENES (Isoprenoides) ---
        "Monoterpene": "CHEBI:25403",      # C10 (Huiles essentielles)
        "Sesquiterpene": "CHEBI:26658",    # C15
        "Diterpene": "CHEBI:23849",        # C20
        "Triterpene": "CHEBI:36326",       # C30
        "Polyprenol": "CHEBI:26197",       # Longues chaînes d'isoprène

        # --- GROUPE 5: MACROCYCLES & CAGES (Topologie unique) ---
        "Porphyrin": "CHEBI:26214",        # Macrocycle 4 pyrroles (Hème)
        "Corrin": "CHEBI:36712",           # Noyau B12 (proche porphyrine)
        "Adamantane": "CHEBI:22268",       # Cage diamant (tricyclique)
        "Crown ether": "CHEBI:46766",      # Éther couronne (cycle O-C-C)
        "Cyclodextrin": "CHEBI:23484",     # Anneau de sucres

        # --- GROUPE 6: GROUPES FONCTIONNELS CHIMIQUES ---
        "Organofluorine": "CHEBI:37143",   # Contient Fluor
        "Organochlorine": "CHEBI:36683",   # Contient Chlore
        "Organobromine": "CHEBI:37141",    # Contient Brome
        "Nitro compound": "CHEBI:35716",   # R-NO2
        "Quinone": "CHEBI:26421",          # Cycle oxydé (C=O conjugué)
        "Epoxide": "CHEBI:32955",          # Cycle à 3 atomes (C-O-C)
        "Lactone": "CHEBI:24973",          # Ester cyclique
        "Lactam": "CHEBI:24971",           # Amide cyclique
        "Nitrile": "CHEBI:18379",          # Groupe cyano -C≡N
        "Isocyanate": "CHEBI:29390",       # Groupe -N=C=O
        "Thiols": "CHEBI:29256"            # Groupe -SH
    }

    print("CREATION LISTE MULTI-CLASSES")

    all_molecules = []
    MOLECULES_PER_CLASS = 100 

    for name, chebi_id in TARGET_CLASSES.items(): 
        target_node = None
        
        if hasattr(ontology, 'get_node'): 
            target_node = ontology.get_node(chebi_id)

        if not target_node:
            continue

        descendants = target_node.get_descendants()
        descendants.add(target_node.chebi_id)

        candidates = list(descendants)
        random.shuffle(candidates) 
        
        count = 0
        for cid in candidates:
            if count >= MOLECULES_PER_CLASS: break
            
            mol_data = db.get_mol(cid)
            if not mol_data: continue
            
            try:
                g = MoleculeGraph.from_moltext(mol_data, chebi_id=cid)
                if len(g.nodes) > 1:
                    all_molecules.append({
                        "chebi_id": cid,
                        "name": f"{name}_{cid}",
                        "graph": g,
                        "class_label": name
                    })
                    count += 1
            except:
                continue
        print(f"Ajouté {count} molécules de {name}.")

    cwl = CWLKernel(similarity="cosine", n_bits=768)
    onto_sim = OntologySimilarity(ontology)
    
    print(f"\nLancement du calcul sur {len(all_molecules)} molécules variées...")
    results = compute_correlation_syntactic_semantic(all_molecules, cwl, onto_sim)
