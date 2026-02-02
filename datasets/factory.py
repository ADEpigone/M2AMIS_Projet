
from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Iterator, Callable, Any, Dict, Union
from dataclasses import dataclass, field
import random
import itertools

import sys
sys.path.append("..")

from datasets.base.dataset import BaseDataset, MoleculeEntry
from datasets.molecule_dataset import MoleculeDataset
from datasets.pair_dataset import PairDataset
from datasets.synthetic_dataset import SyntheticDatasetConfig, SyntheticMoleculeGenerator
from datasets.dataloader import DataLoader

class DatasetFactory:
    """
    Factory pour créer des datasets facilement.
    """
    
    @staticmethod
    def from_chebi_ids(
        chebi_ids: List[str],
        name: str = "chebi_dataset",
        chebi_retriever=None
    ) -> MoleculeDataset:
        """
        Crée un dataset à partir d'une liste de ChEBI IDs.

        """
        molecules = []
        
        for chebi_id in chebi_ids:
            try:
                if chebi_retriever is not None:
                    from graph import MoleculeGraph
                    mol_text = chebi_retriever.get_mol(chebi_id)
                    graph = MoleculeGraph.from_moltext(mol_text)
                else:
                    graph = None
                
                mol_entry = MoleculeEntry(
                    chebi_id=chebi_id,
                    graph=graph,
                    properties={"source": "chebi"}
                )
                molecules.append(mol_entry)
            except Exception as e:
                print(f"Warning: Could not load molecule {chebi_id}: {e}")
        
        return MoleculeDataset(molecules, name=name)
    
    @staticmethod
    def from_mol_files(
        file_paths: List[str],
        name: str = "local_dataset"
    ) -> MoleculeDataset:
        """
        Crée un dataset à partir de fichiers mol locaux.
        """
        import os
        molecules = []
        
        for path in file_paths:
            try:
                from graph import MoleculeGraph
                graph = MoleculeGraph.from_molfile(path)
                filename = os.path.basename(path)
                mol_id = filename.replace(".mol", "")
                
                mol_entry = MoleculeEntry(
                    chebi_id=mol_id,
                    graph=graph,
                    properties={"source": "local", "file_path": path}
                )
                molecules.append(mol_entry)
            except Exception as e:
                print(f"Warning: Could not load file {path}: {e}")
        
        return MoleculeDataset(molecules, name=name)
    
    @staticmethod
    def create_synthetic(
        config: Optional[SyntheticDatasetConfig] = None,
        **kwargs
    ) -> Tuple[MoleculeDataset, PairDataset]:
        """
        Crée un dataset synthétique.
        """
        if config is None:
            config = SyntheticDatasetConfig(**kwargs)
        
        generator = SyntheticMoleculeGenerator(config)
        return generator.generate_dataset()
    
    @staticmethod
    def create_benchmark_dataset(
        num_small: int = 50,
        num_medium: int = 30,
        num_large: int = 20,
        seed: int = 42
    ) -> MoleculeDataset:
        """
        Crée un dataset de benchmark avec différentes tailles de molécules.
        """
        all_molecules = []
        
        # Petites molécules
        small_config = SyntheticDatasetConfig(
            num_molecules=num_small,
            min_atoms=3,
            max_atoms=10,
            seed=seed,
            include_isomorphic_pairs=False
        )
        small_gen = SyntheticMoleculeGenerator(small_config)
        small_ds, _ = small_gen.generate_dataset()
        for mol in small_ds:
            mol.properties["size_category"] = "small"
            all_molecules.append(mol)
        
        # Molécules moyennes
        medium_config = SyntheticDatasetConfig(
            num_molecules=num_medium,
            min_atoms=10,
            max_atoms=30,
            seed=seed + 1,
            include_isomorphic_pairs=False
        )
        medium_gen = SyntheticMoleculeGenerator(medium_config)
        medium_ds, _ = medium_gen.generate_dataset()
        for mol in medium_ds:
            mol.properties["size_category"] = "medium"
            mol.chebi_id = mol.chebi_id.replace("SYNTHETIC_", "SYNTHETIC_MED_")
            all_molecules.append(mol)
        
        # Grandes molécules
        large_config = SyntheticDatasetConfig(
            num_molecules=num_large,
            min_atoms=30,
            max_atoms=50,
            seed=seed + 2,
            include_isomorphic_pairs=False
        )
        large_gen = SyntheticMoleculeGenerator(large_config)
        large_ds, _ = large_gen.generate_dataset()
        for mol in large_ds:
            mol.properties["size_category"] = "large"
            mol.chebi_id = mol.chebi_id.replace("SYNTHETIC_", "SYNTHETIC_LRG_")
            all_molecules.append(mol)
        
        return MoleculeDataset(all_molecules, name="benchmark_dataset")
    

if __name__ == "__main__":
    print("Datasets Module Demo")
    # Créer un dataset synthétique
    config = SyntheticDatasetConfig(
        num_molecules=15,
        min_atoms=5,
        max_atoms=12,
        seed=42,
        include_isomorphic_pairs=True,
        num_isomorphic_pairs=3
    )
    mol_dataset, pair_dataset = DatasetFactory.create_synthetic(config)
    print(f"\nDataset créé: {mol_dataset}")
    print(f"Paires ground truth: {len(pair_dataset)}")
    # DataLoader
    loader = DataLoader(mol_dataset, batch_size=4, shuffle=True)
    print(f"\nDataLoader avec batch_size=4: {len(loader)} batches")
    # Afficher quelques molécules
    print("\nPremières molécules:")
    for i, mol in enumerate(mol_dataset):
        if i >= 3:
            break
        print(f"  {mol}")
    # Afficher les paires isomorphes
    print("\nPaires isomorphes (ground truth):")
    for pair in pair_dataset:
        print(f"  {pair.mol1.chebi_id} <-> {pair.mol2.chebi_id}")
    print("\n" + "=" * 60)
    print("Demo terminée!")
