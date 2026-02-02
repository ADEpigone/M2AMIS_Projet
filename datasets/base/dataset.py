from abc import ABC, abstractmethod
from typing import Iterator, Any, Dict, Optional
from dataclasses import dataclass, field
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
