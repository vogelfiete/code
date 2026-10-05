"""The data types the xlms package hands back: a single crosslink and a
whole dataset of crosslinks plus protein sequences."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


# -----------------------------------------------------------------------------
# One crosslink
#
# Two proteins, the linked position in each, the score and the target/decoy
# classification of the link.
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# A whole dataset
#
# All crosslinks of a file together with the protein sequences, plus a few
# shortcuts to pick out the links you usually want.
# -----------------------------------------------------------------------------

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
