from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Iterator, Callable, Any, Dict, Union
from dataclasses import dataclass, field
import random
import itertools

from datasets.base.dataset import BaseDataset


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