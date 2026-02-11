import numpy as np
from clustering.clusters import DB_PATH, load_clusters
from utils import load_ontology
from Chebi.CheBi2 import CheBi2
def analyze_family_purity(family_name_or_id: str, ontology, molecules: list, cluster_map: dict):
    """
    Analyse la pureté des clusters pour une famille donnée de l'ontologie.
    
    Args:
        family_name_or_id: Nom (ex: 'peptide') ou ID Chebi (ex: 'CHEBI:16670')
        ontology: L'objet ontologie chargé
        molecules: La liste de toutes les molécules chargées (avec chebi_id)
        cluster_map: Dict {cluster_id: [liste_molecules]}
    """
    print(f"ANALYSE DE LA FAMILLE : {family_name_or_id.upper()}")

    root_node = None
    if family_name_or_id.startswith("CHEBI:"):
        root_node = ontology.get_node(family_name_or_id)
    else:
        # Recherche par nom insensible à la casse
        target = family_name_or_id.lower()
        for node in ontology.nodes.values():
            if node.name and node.name.lower() == target:
                root_node = node
                break
    
    if not root_node:
        print(f"Famille '{family_name_or_id}' introuvable dans l'ontologie.")
        return

    print(f"Racine trouvée: {root_node.name} ({root_node.chebi_id})")

    family_ids = {root_node.chebi_id}
    descendants = root_node.get_descendants() 
    if descendants:
        family_ids.update(descendants)
    
    print(f"Total de concepts dans la famille: {len(family_ids)}")

    target_mols_ids = set()
    target_mols_count = 0
    
    for mol in molecules:
        cid = mol[0]
        
        mol_node = ontology.get_node(cid)
        if not mol_node: continue

        ancestors = mol_node.get_ancestors()
        if (cid in family_ids) or (ancestors & family_ids):
            target_mols_ids.add(str(mol[0])) # On stocke en string pour comparer facilement
            target_mols_count += 1

    pct_global = 100 * target_mols_count / len(molecules) if molecules else 0
    print(f"Molécules trouvées: {target_mols_count}/{len(molecules)} ({pct_global:.1f}%)")

    if target_mols_count == 0:
        print("Aucune molécule de cette famille dans le dataset.")
        return

    cluster_stats = [] # (cluster_id, total_size, count_family, ratio)
    
    for cid, members in cluster_map.items():
        #print(members[0])
        count = sum(1 for m in members if str(m["chebi_id"]) in target_mols_ids)
        if count > 0:
            ratio = count / len(members)
            cluster_stats.append((cid, len(members), count, ratio))
    
    cluster_stats.sort(key=lambda x: x[3], reverse=True)

    pure_clusters = [c for c in cluster_stats if c[3] == 1.0 and c[1] > 1]
    mixed_clusters = [c for c in cluster_stats if 0.5 <= c[3] < 1.0]

    print(f"\nDISTRIBUTION DES CLUSTERS ({len(cluster_stats)} concernés)")
    print(f"\t- Clusters 100% purs (>1 élément) : {len(pure_clusters)}")
    print(f"\t- Clusters Mixtes (50-99%)        : {len(mixed_clusters)}")

    if pure_clusters:
        print(f"\nTop 5 Clusters PURS (100% {family_name_or_id}):")
        pure_clusters.sort(key=lambda x: x[1], reverse=True)
        for cid, size, count, ratio in pure_clusters[:5]:
            members = cluster_map[cid]
            names = [m["name"] for m in members[:3]]
            print(f"\t• Cluster {cid} (n={size}): {', '.join(names)}...")

    if mixed_clusters:
        print(f"\nTop 5 Clusters MIXTES (pollués):")
        mixed_clusters.sort(key=lambda x: x[1], reverse=True)
        for cid, size, count, ratio in mixed_clusters[:5]:
            print(f"\t• Cluster {cid} (n={size}): {count} cibles ({ratio:.1%})")

    ratios = [c[3] for c in cluster_stats]
    print(f"\nMétriques de Pureté (sur les clusters contenant au moins 1 cible):")
    print(f"\t- Moyenne : {np.mean(ratios):.1%}")
    print(f"\t- Médiane : {np.median(ratios):.1%}")

if __name__ == "__main__":
    cluster_map = load_clusters()
    ontology = load_ontology()

    chebi_db = CheBi2(DB_PATH)
    molecules = list(chebi_db.get_all_mols())

    analyze_family_purity("flavonoid", ontology, molecules, cluster_map)