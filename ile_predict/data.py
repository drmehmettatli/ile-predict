"""Structure input: resolve names to SMILES (PubChem) and canonicalize (RDKit)."""
from __future__ import annotations
from typing import Optional
import time


def name_to_smiles(name: str, retries: int = 3, pause: float = 2.0) -> Optional[str]:
    """Resolve a compound name to an isomeric SMILES via PubChem. None on failure."""
    import pubchempy as pcp
    for _ in range(retries):
        try:
            hits = pcp.get_compounds(name, "name")
            return hits[0].smiles if hits else None
        except Exception:
            time.sleep(pause)
    return None


def canonicalize(smiles: str) -> Optional[str]:
    """Return RDKit canonical SMILES, or None if unparseable."""
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def applicability_flag(mw: float) -> str:
    """Flag molecules outside ADMET-AI's drug-like training domain."""
    if mw is None:
        return "unknown"
    if mw > 600 or mw < 100:
        return "domain-edge (prediction unreliable)"
    return "ok"
