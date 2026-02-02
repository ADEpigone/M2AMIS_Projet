from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Iterator, Callable, Any, Dict, Union
from dataclasses import dataclass, field
import random
import itertools

from datasets.base.dataset import BaseDataset, MoleculeEntry

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