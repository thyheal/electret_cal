from rdkit import Chem
from rdkit.Chem import AllChem
def FFKM_builder(m,n,c,connection_direction):
    FFKM = n * 'C(F)(F)' +  m * 'C(F)(F)C(OC(F)(F)(F))(F)'
    if connection_direction == 'L':
        connection_point = 0
        FFKM = c*'C'+ FFKM + 'F'
    elif connection_direction == 'R':
        FFKM = 'F' + FFKM + 'C' * c
        if c == 0:
            connection_point = 3 * n + 10 * (m - 1) + 3 + 1
        else:
            connection_point = 3 * n + 10 * m + c
    return FFKM, connection_point


def crosslink_reaction(base,FFKM,connection_point):
    mol = Chem.MolFromSmiles(base)
    patt = Chem.MolFromSmarts('[C;!R]C=[C;!R]')
    repl = Chem.MolFromSmiles(FFKM)
    rms = AllChem.ReplaceSubstructs(mol, patt, repl, replaceAll=True, replacementConnectionPoint=connection_point)
    s1 = Chem.MolToSmiles(rms[0])
    return s1