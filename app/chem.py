from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from app.config import TEMP_DIR


def process_smiles(smi: str):
    """Read smiles and return canon_smiles, canon_mb, temp_path"""
    m = Chem.MolFromSmiles(smi)
    fpath = TEMP_DIR.joinpath("temp.mol")
    m_canon = rdkit_atom_order(m)
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


def canonicalize_atom_order(m, reverse=True, add_hs=True):
    """Canonicalize using RDKIT

    Args:
        m (rdkit.Chem.Mol): Mol object for RDKit

    Returns:
        rdkit.Chem.Mol: New canonicalized RDKit mol
    """
    if add_hs:
        mH = Chem.AddHs(m)
    else:
        mH = m
    Compute2DCoords(mH)
    m_neworder = tuple(
        zip(
            *sorted(
                [(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mH))],
                reverse=reverse,
            )
        )
    )[1]
    m_canon = Chem.RenumberAtoms(mH, m_neworder)
    add_atom_indices(m_canon)
    return m_canon


def rdkit_atom_order(m, add_hs=True):
    """Canonicalize using RDKit SMILES export

    Args:
        m (rdkit.Chem.Mol): Mol object for RDKit

    Returns:
        rdkit.Chem.Mol: New canonicalized RDKit mol
    """
    m_renum = Chem.MolFromSmiles(Chem.MolToSmiles(m))
    if add_hs:
        m_canon = Chem.AddHs(m_renum)
    else:
        m_canon = m_renum
    add_atom_indices(m_canon)
    return m_canon


def add_atom_indices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i + 1)


def smi_to_mol(smi):
    return Chem.MolFromSmiles(smi)


def mol_to_smi(m):
    return Chem.MolToSmiles(m)
