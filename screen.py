from rdkit import Chem
from rdkit.Chem import AllChem

def filter(smiles,functional_group):
    mol = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts(functional_group)
    matches = mol.GetSubstructMatches(patt)
    if matches:
        return True
    else:
        return False