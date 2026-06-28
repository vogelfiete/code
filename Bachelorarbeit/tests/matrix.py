import os
import sys
import tempfile
import textwrap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from xlms import CrossLink, CrossLinkDataset, load_dataset, read_csv, read_fasta

_SAMPLE_CSV = textwrap.dedent("""\
    Id,Protein1,Protein2,SeqPos1,SeqPos2,Score,isDecoy,isTT,isTD,isDD
    1,PROT_A,PROT_B,10,25,42.7,False,True,False,False
    2,PROT_A,DECOY_B,10,25,18.3,True,False,True,False
    3,DECOY_A,DECOY_B,5,12,9.1,True,False,False,True
    4,PROT_A,PROT_B,7,30,55.0,0,1,0,0
""")

_SAMPLE_FASTA = textwrap.dedent("""\
    >PROT_A
    MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK
    >PROT_B
    MSTAGKVIKCKAAVLWEEKKPFSIEEVEVAPPKAHEVRIKMVATGICRSDDHVV
    >DECOY_A
    LLSPSVAERFSVAAGLCAFSMPIKPDLQEVLAGICATMVVMFCNIKAGLREDGV
    >DECOY_B
    QFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDA
""")


def _write_temp(content: str, suffix: str) -> str:
    f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False, encoding="utf-8")
    f.write(content)
    f.close()
    return f.name


def test_read_csv_types():
    path = _write_temp(_SAMPLE_CSV, ".csv")
    links = read_csv(path)
    os.unlink(path)

    assert len(links) == 4, f"Expected 4 crosslinks, got {len(links)}"
    first = links[0]
    assert isinstance(first, CrossLink)
    assert first.id == "1"
    assert first.protein1 == "PROT_A"
    assert first.protein2 == "PROT_B"
    assert first.seq_pos1 == 10
    assert first.seq_pos2 == 25
    assert isinstance(first.score, float)
    assert type(first.is_decoy) is bool
    assert type(first.is_tt) is bool
    assert first.is_decoy is False
    assert first.is_tt is True

    second = links[1]
    assert second.is_decoy is True
    assert second.is_td is True

    # Row 4 uses 0/1 integers — verify they coerce to bool
    fourth = links[3]
    assert fourth.is_decoy is False
    assert fourth.is_tt is True
    print("test_read_csv_types PASSED")


def test_read_fasta_types():
    path = _write_temp(_SAMPLE_FASTA, ".fasta")
    seqs = read_fasta(path)
    os.unlink(path)

    assert set(seqs.keys()) == {"PROT_A", "PROT_B", "DECOY_A", "DECOY_B"}
    assert isinstance(seqs["PROT_A"], str)
    assert seqs["PROT_A"][0] == "M"
    assert seqs["PROT_A"] == seqs["PROT_A"].upper()
    print("test_read_fasta_types PASSED")


def test_load_dataset_filters():
    csv_path = _write_temp(_SAMPLE_CSV, ".csv")
    fasta_path = _write_temp(_SAMPLE_FASTA, ".fasta")
    ds = load_dataset(csv_path, fasta_path)
    os.unlink(csv_path)
    os.unlink(fasta_path)

    assert isinstance(ds, CrossLinkDataset)
    assert len(ds.crosslinks) == 4
    assert len(ds.sequences) == 4
    assert len(ds.target_target()) == 2
    assert len(ds.target_decoy()) == 1
    assert len(ds.decoy_decoy()) == 1
    assert len(ds.above_score(40.0)) == 2
    print("test_load_dataset_filters PASSED")


_ALTERNATE_CSV = textwrap.dedent("""\
    Id,Protein1,Protein2,SeqPos1,SeqPos2,Score,Decoy1,Decoy2,DecoyType
    1,PROT_A,PROT_B,10,25,42.7,False,False,TT
    2,PROT_A,DECOY_B,10,25,18.3,False,True,TD
    3,DECOY_A,DECOY_B,5,12,9.1,True,True,DD
    4,PROT_A,PROT_B,7,30,55.0,0,0,TT
""")

_BROKEN_CSV = textwrap.dedent("""\
    Id,Protein1,Protein2,SeqPos1,SeqPos2,Score,SomeOtherCol
    1,PROT_A,PROT_B,10,25,42.7,xyz
""")


def test_read_csv_alternate_format():
    path = _write_temp(_ALTERNATE_CSV, ".csv")
    links = read_csv(path)
    os.unlink(path)

    assert len(links) == 4
    first = links[0]
    assert first.is_decoy is False
    assert first.is_tt is True
    assert first.is_td is False
    assert first.is_dd is False

    second = links[1]
    assert second.is_decoy is True   # Decoy2=True
    assert second.is_td is True

    third = links[2]
    assert third.is_decoy is True
    assert third.is_dd is True

    # Row 4 uses 0/0 integers for Decoy1/Decoy2
    fourth = links[3]
    assert fourth.is_decoy is False
    assert fourth.is_tt is True

    # Missing required columns should raise
    broken_path = _write_temp(_BROKEN_CSV, ".csv")
    try:
        read_csv(broken_path)
        assert False, "Expected ValueError"
    except ValueError:
        pass
    finally:
        os.unlink(broken_path)

    print("test_read_csv_alternate_format PASSED")


if __name__ == "__main__":
    test_read_csv_types()
    test_read_fasta_types()
    test_load_dataset_filters()
    test_read_csv_alternate_format()
    print("\nAll smoke tests passed.")
