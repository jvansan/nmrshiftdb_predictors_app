from rdkit import Chem

from app.config import TEMP_DIR


def process_smiles(smi: str):
    """Read smiles and return canon_smiles, canon_mb, temp_path"""
    m = Chem.MolFromSmiles(smi)
    fpath = TEMP_DIR.joinpath("temp.mol")
    m_canon = canonicalize_atom_order(m)
    if not TEMP_DIR.is_dir():
        TEMP_DIR.mkdir()
    Chem.MolToMolFile(m_canon, str(fpath))
    return (Chem.MolToSmiles(m_canon), Chem.MolToMolBlock(m_canon), fpath)


def process_inchi(inchi: str):
    """Read inchi and return canon_smiles, canon_mb, temp_path"""
    m = Chem.MolFromInchi(inchi)
    fpath = TEMP_DIR.joinpath("temp.mol")
    m_canon = canonicalize_atom_order(m)
    if not TEMP_DIR.is_dir():
        TEMP_DIR.mkdir()
    Chem.MolToMolFile(m_canon, str(fpath))
    return (Chem.MolToSmiles(m_canon), Chem.MolToMolBlock(m_canon), fpath)


def canonicalize_atom_order(m):
    m = Chem.AddHs(m)
    m_neworder = tuple(
        zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(m))]))
    )[1]
    m_renum = Chem.RenumberAtoms(m, m_neworder)
    add_atom_indices(m_renum)
    return m_renum


def add_atom_indices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i + 1)
