
from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Iterator, Callable, Any, Dict, Union
from dataclasses import dataclass, field
import random
import itertools

from datasets.base.dataset import BaseDataset, MoleculeEntry
from datasets.molecule_dataset import MoleculeDataset

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