
from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Iterator, Callable, Any, Dict, Union
from dataclasses import dataclass, field
import random
import itertools

from datasets.base.dataset import BaseDataset, MoleculeEntry
from datasets.molecule_dataset import MoleculeDataset
from datasets.pair_dataset import PairDataset, MoleculePair

@dataclass
class SyntheticDatasetConfig:
    """
    Configuration pour la génération de datasets synthétiques.
    
    Attributes:
        num_molecules: Nombre de molécules à générer
        min_atoms: Nombre minimum d'atomes
        max_atoms: Nombre maximum d'atomes
        atom_types: Types d'atomes possibles
        bond_types: Types de liaisons possibles
        connectivity: Degré de connectivité (0.0 à 1.0)
        seed: Graine aléatoire pour reproductibilité
        include_isomorphic_pairs: Générer des paires isomorphes connues
        num_isomorphic_pairs: Nombre de paires isomorphes à générer
    """
    num_molecules: int = 100
    min_atoms: int = 3
    max_atoms: int = 20
    atom_types: List[str] = field(default_factory=lambda: ["C", "H", "O", "N", "S", "P"])
    bond_types: List[str] = field(default_factory=lambda: ["1", "2", "3"])
    connectivity: float = 0.3
    seed: Optional[int] = None
    include_isomorphic_pairs: bool = True
    num_isomorphic_pairs: int = 10


class SyntheticMoleculeGenerator:
    """
    Générateur de molécules synthétiques.

    """
    
    def __init__(self, config: SyntheticDatasetConfig):
        self.config = config
        if config.seed is not None:
            random.seed(config.seed)
    
    def _generate_single_molecule(self, mol_id: str) -> MoleculeEntry:
        """
        Génère une molécule synthétique aléatoire.
        """
        num_atoms = random.randint(self.config.min_atoms, self.config.max_atoms)
        
        # Générer les noeuds (atomes)
        nodes = []
        for i in range(num_atoms):
            atom_type = random.choice(self.config.atom_types)
            nodes.append({"id": i + 1, "color": atom_type})
        
        # Générer les arêtes (liaisons)
        edges = []
        # D'abord, créer un arbre couvrant pour assurer la connectivité
        for i in range(1, num_atoms):
            j = random.randint(0, i - 1)
            bond_type = random.choice(self.config.bond_types)
            edges.append({"u": j + 1, "v": i + 1, "color": bond_type})
        
        # Ajouter des arêtes supplémentaires selon la connectivité
        max_additional = int(num_atoms * (num_atoms - 1) / 2 - (num_atoms - 1))
        num_additional = int(max_additional * self.config.connectivity)
        
        existing_edges = set((min(e["u"], e["v"]), max(e["u"], e["v"])) for e in edges)
        
        for _ in range(num_additional):
            attempts = 0
            while attempts < 10:
                i = random.randint(1, num_atoms)
                j = random.randint(1, num_atoms)
                if i != j and (min(i, j), max(i, j)) not in existing_edges:
                    bond_type = random.choice(self.config.bond_types)
                    edges.append({"u": i, "v": j, "color": bond_type})
                    existing_edges.add((min(i, j), max(i, j)))
                    break
                attempts += 1
        
        graph_data = {"nodes": nodes, "edges": edges}
        
        return MoleculeEntry(
            chebi_id=f"SYNTHETIC_{mol_id}",
            graph=graph_data,
            name=f"Synthetic Molecule {mol_id}",
            properties={
                "num_atoms": num_atoms,
                "num_bonds": len(edges),
                "is_synthetic": True
            }
        )
    
    def _create_isomorphic_copy(self, original: MoleculeEntry, new_id: str) -> MoleculeEntry:
        """
        Crée une copie isomorphe d'une molécule (permutation des noeuds).
        """
        graph_data = original.graph
        nodes = graph_data["nodes"]
        edges = graph_data["edges"]
        
        # Créer une permutation aléatoire des IDs
        old_ids = [n["id"] for n in nodes]
        new_ids = old_ids.copy()
        random.shuffle(new_ids)
        id_mapping = dict(zip(old_ids, new_ids))
        
        # Appliquer la permutation
        new_nodes = [{"id": id_mapping[n["id"]], "color": n["color"]} for n in nodes]
        new_edges = [{"u": id_mapping[e["u"]], "v": id_mapping[e["v"]], "color": e["color"]} for e in edges]
        
        return MoleculeEntry(
            chebi_id=f"SYNTHETIC_{new_id}",
            graph={"nodes": new_nodes, "edges": new_edges},
            name=f"Synthetic Molecule {new_id} (isomorphic to {original.chebi_id})",
            properties={
                **original.properties,
                "isomorphic_to": original.chebi_id
            }
        )
    
    def generate_dataset(self) -> Tuple[MoleculeDataset, PairDataset]:
        """
        Génère un dataset synthétique complet.

        """
        molecules = []
        isomorphic_pairs = []
        
        # Générer les molécules de base
        base_count = self.config.num_molecules
        if self.config.include_isomorphic_pairs:
            base_count -= self.config.num_isomorphic_pairs
        
        for i in range(base_count):
            mol = self._generate_single_molecule(str(i))
            molecules.append(mol)
        
        # Générer des paires isomorphes
        if self.config.include_isomorphic_pairs:
            for i in range(self.config.num_isomorphic_pairs):
                original = random.choice(molecules[:base_count])
                iso_copy = self._create_isomorphic_copy(original, f"ISO_{i}")
                molecules.append(iso_copy)
                isomorphic_pairs.append(MoleculePair(
                    mol1=original,
                    mol2=iso_copy,
                    is_isomorphic=True
                ))
        
        mol_dataset = MoleculeDataset(molecules, name="synthetic_dataset")
        pair_dataset = PairDataset(isomorphic_pairs, name="synthetic_pairs_ground_truth")
        
        return mol_dataset, pair_dataset