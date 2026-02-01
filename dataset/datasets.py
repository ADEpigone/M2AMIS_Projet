"""
Datasets Module
Ce module contient les classes Dataset et DataLoader pour le projet.

Classes:
    - BaseDataset: Classe abstraite de base
    - MoleculeEntry: Dataclass pour une molécule
    - MoleculeDataset: Dataset de molécules
    - MoleculePair: Dataclass pour une paire de molécules
    - PairDataset: Dataset de paires
    - SyntheticDatasetConfig: Configuration pour génération synthétique
    - SyntheticMoleculeGenerator: Générateur de molécules synthétiques
    - DataLoader: Pour itérer avec batching
    - DatasetFactory: Factory pour créer des datasets facilement
"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Iterator, Callable, Any, Dict, Union
from dataclasses import dataclass, field
import random
import itertools

# BASE CLASSES

class BaseDataset(ABC):
    """
    Classe abstraite de base pour tous les datasets.
    Inspirée de torch.utils.data.Dataset.
    
    Tout dataset doit implémenter:
    - __len__: retourne la taille du dataset
    - __getitem__: retourne un élément par index
    """
    
    @abstractmethod
    def __len__(self) -> int:
        pass
    
    @abstractmethod
    def __getitem__(self, index: int):
        pass
    
    def __iter__(self) -> Iterator:
        for i in range(len(self)):
            yield self[i]
    
    def get_metadata(self) -> Dict[str, Any]:
        return {
            "type": self.__class__.__name__,
            "size": len(self)
        }

# MOLECULE ENTRY

@dataclass
class MoleculeEntry:
    """
    Entrée représentant une molécule dans le dataset.
    Attributes:
        chebi_id: Identifiant ChEBI de la molécule
        graph: Le graphe de la molécule (MoleculeGraph)
        name: Nom de la molécule (optionnel)
        properties: Propriétés additionnelles
    """
    chebi_id: str
    graph: Any
    name: Optional[str] = None
    properties: Dict[str, Any] = field(default_factory=dict)
    
    def __hash__(self):
        return hash(self.chebi_id)
    
    def __eq__(self, other):
        if isinstance(other, MoleculeEntry):
            return self.chebi_id == other.chebi_id
        return False
    
    def __repr__(self):
        return f"MoleculeEntry(chebi_id='{self.chebi_id}', name='{self.name}')"


# MOLECULE DATASET

class MoleculeDataset(BaseDataset):
    """
    Dataset contenant des molécules.
    
    Peut être créé:
    - À partir d'une liste de ChEBI IDs
    - À partir de fichiers mol locaux
    - À partir d'un filtrage par propriétés
    """
    
    def __init__(
        self,
        molecules: Optional[List[MoleculeEntry]] = None,
        name: str = "unnamed_dataset"
    ):
        """
        Initialise le dataset.
        
        Args:
            molecules: Liste des entrées de molécules
            name: Nom du dataset
        """
        self._molecules: List[MoleculeEntry] = molecules or []
        self._name = name
        self._index_map: Dict[str, int] = {}
        self._rebuild_index()
    
    def _rebuild_index(self):
        self._index_map = {mol.chebi_id: i for i, mol in enumerate(self._molecules)}
    
    def __len__(self) -> int:
        return len(self._molecules)
    
    def __getitem__(self, index: Union[int, str]) -> MoleculeEntry:
        #Récupère une molécule par index ou par chebi_id.

        if isinstance(index, str):
            if index not in self._index_map:
                raise KeyError(f"Molecule with chebi_id '{index}' not found")
            return self._molecules[self._index_map[index]]
        return self._molecules[index]
    
    def add_molecule(self, molecule: MoleculeEntry):
        """Ajoute une molécule au dataset."""
        if molecule.chebi_id in self._index_map:
            raise ValueError(f"Molecule {molecule.chebi_id} already exists")
        self._index_map[molecule.chebi_id] = len(self._molecules)
        self._molecules.append(molecule)
    
    def remove_molecule(self, chebi_id: str):
        """Retire une molécule du dataset."""
        if chebi_id not in self._index_map:
            raise KeyError(f"Molecule {chebi_id} not found")
        idx = self._index_map[chebi_id]
        self._molecules.pop(idx)
        self._rebuild_index()
    
    def filter(self, predicate: Callable[[MoleculeEntry], bool]) -> 'MoleculeDataset':
        """
        Filtre le dataset selon un prédicat.

        """
        filtered = [mol for mol in self._molecules if predicate(mol)]
        return MoleculeDataset(filtered, name=f"{self._name}_filtered")
    
    def get_all_chebi_ids(self) -> List[str]:
        #Retourne tous les ChEBI IDs.
        return list(self._index_map.keys())
    
    def get_pairs(self, include_self: bool = False) -> Iterator[Tuple[MoleculeEntry, MoleculeEntry]]:
        if include_self:
            return itertools.combinations_with_replacement(self._molecules, 2)
        return itertools.combinations(self._molecules, 2)
    
    def get_metadata(self) -> Dict[str, Any]:
        # Retourne les métadonnées du dataset.
        base = super().get_metadata()
        base.update({
            "name": self._name,
            "chebi_ids": self.get_all_chebi_ids()
        })
        return base
    
    def __repr__(self):
        return f"MoleculeDataset(name='{self._name}', size={len(self)})"


# PAIR DATASET

@dataclass
class MoleculePair:
    """
    Représente une paire de molécules avec un label optionnel.
    Attributes:
        mol1: Première molécule
        mol2: Deuxième molécule
        is_isomorphic: Label indiquant si les molécules sont isomorphes
        similarity: Score de similarité (optionnel)
    """
    mol1: MoleculeEntry
    mol2: MoleculeEntry
    is_isomorphic: Optional[bool] = None
    similarity: Optional[float] = None


class PairDataset(BaseDataset):
    """
    Dataset de paires de molécules.
    Utile pour le benchmarking des algorithmes d'isomorphisme.
    """
    
    def __init__(
        self,
        pairs: Optional[List[MoleculePair]] = None,
        name: str = "unnamed_pair_dataset"
    ):
        self._pairs = pairs or []
        self._name = name
    
    def __len__(self) -> int:
        return len(self._pairs)
    
    def __getitem__(self, index: int) -> MoleculePair:
        return self._pairs[index]
    
    def add_pair(self, pair: MoleculePair):
        self._pairs.append(pair)
    
    @classmethod
    def from_molecule_dataset(
        cls,
        dataset: MoleculeDataset,
        include_self: bool = False
    ) -> 'PairDataset':
        pairs = [
            MoleculePair(mol1, mol2)
            for mol1, mol2 in dataset.get_pairs(include_self)
        ]
        return cls(pairs, name=f"{dataset._name}_pairs")
    
    def get_positive_pairs(self) -> List[MoleculePair]:
        """Retourne les paires isomorphes (ground truth)."""
        return [p for p in self._pairs if p.is_isomorphic is True]
    
    def get_negative_pairs(self) -> List[MoleculePair]:
        """Retourne les paires non-isomorphes."""
        return [p for p in self._pairs if p.is_isomorphic is False]
    
    def get_metadata(self) -> Dict[str, Any]:
        base = super().get_metadata()
        base.update({
            "name": self._name,
            "num_positive": len(self.get_positive_pairs()),
            "num_negative": len(self.get_negative_pairs())
        })
        return base
    
    def __repr__(self):
        return f"PairDataset(name='{self._name}', size={len(self)})"

# SYNTHETIC DATASET GENERATION

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


# DATALOADER

class DataLoader:
    """
    DataLoader pour itérer sur les datasets.
    Supporte le batching et le shuffling
    """
    
    def __init__(
        self,
        dataset: BaseDataset,
        batch_size: int = 1,
        shuffle: bool = False,
        drop_last: bool = False
    ):
        self.dataset = dataset
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.drop_last = drop_last
    
    def __len__(self) -> int:
        """Nombre de batches."""
        if self.drop_last:
            return len(self.dataset) // self.batch_size
        return (len(self.dataset) + self.batch_size - 1) // self.batch_size
    
    def __iter__(self) -> Iterator[List]:
        """Itère sur les batches."""
        indices = list(range(len(self.dataset)))
        
        if self.shuffle:
            random.shuffle(indices)
        
        batch = []
        for idx in indices:
            batch.append(self.dataset[idx])
            if len(batch) == self.batch_size:
                yield batch
                batch = []
        
        if batch and not self.drop_last:
            yield batch

# DATASET FACTORY

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


# EXPORTS

__all__ = [
    'BaseDataset',
    'MoleculeEntry',
    'MoleculeDataset',
    'MoleculePair',
    'PairDataset',
    'SyntheticDatasetConfig',
    'SyntheticMoleculeGenerator',
    'DataLoader',
    'DatasetFactory'
]

# DEMO

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
