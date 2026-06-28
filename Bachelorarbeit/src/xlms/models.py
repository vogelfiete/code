from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


@dataclass(frozen=True, slots=True)
class CrossLink:
    id: str
    protein1: str
    protein2: str
    seq_pos1: int
    seq_pos2: int
    score: float
    is_decoy: bool
    is_tt: bool
    is_td: bool
    is_dd: bool


@dataclass
class CrossLinkDataset:
    crosslinks: List[CrossLink]
    sequences: Dict[str, str]

    def target_target(self) -> List[CrossLink]:
        return [xl for xl in self.crosslinks if xl.is_tt]

    def target_decoy(self) -> List[CrossLink]:
        return [xl for xl in self.crosslinks if xl.is_td]

    def decoy_decoy(self) -> List[CrossLink]:
        return [xl for xl in self.crosslinks if xl.is_dd]

    def above_score(self, threshold: float) -> List[CrossLink]:
        return [xl for xl in self.crosslinks if xl.score >= threshold]
