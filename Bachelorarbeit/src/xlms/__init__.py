"""Loading XL-MS crosslink data: CSV files with crosslinks and FASTA files
with the protein sequences."""
from xlms.io import load_dataset, read_csv, read_fasta
from xlms.models import CrossLink, CrossLinkDataset

__all__ = [
    "CrossLink",
    "CrossLinkDataset",
    "read_csv",
    "read_fasta",
    "load_dataset",
]
