# import utils
# import gaussian_cal
# import pandas as pd
# import subprocess
# import os
# import time
# from rdkit.Chem import PandasTools
# import pandas as pd
# from tqdm import tqdm
# import seaborn as sns
# import numpy as np
# import matplotlib.pyplot as plt
# import scipy.stats as stats
# from rdkit import Chem
# from rdkit.Chem import AllChem

print(1)

name_list = ['A13931', 'A94974', 'A15322697', 'A12612564', 'A15322699', 'A15051391', 'A12612565', 'A167189779', 'A163913855', 'A163511708', 'A157680575', 'A155767072', 'A152868161', 'A146274067', 'A145235579', 'A144512334', 'A140011603', 'A138624894', 'A137495078', 'A134296113', 'A132041350', 'A118198862', 'A90417328', 'A90337160', 'A90265042', 'A90243281', 'A89902120', 'A89716277', 'A89562383', 'A89305370', 'A89269359', 'A89169074', 'A89141817', 'A68747439', 'A59116104', 'A58472088', 'A54347993', 'A45844598', 'A23042185', 'A23042184', 'A20188883', 'A15620668', 'A15414735', 'A12455325', 'A60341562', 'A101375506', 'A130316971', 'A137495079']
smiles_list = ['O=c1n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'O=c1[nH]c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'Cn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'Cn1c(=O)n(C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=C1N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)CN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'CCCCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=c1n(O)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'CC(C)n1c(=O)n(C(C)C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCCCN1C(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)CN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C1=O', 'C=Cn1c(=O)n(CC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCCCn1c(=O)n(CCCC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'C#CCn1c(=O)n(CC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'C=C1N(C)C(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'CCCn1c(=O)n(CC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCCCn1c(=O)n(C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=c1n(CC(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'COCn1c(=O)n(COC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=C1N(F)CN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'Cn1c(=O)n(CC(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CC(C)n1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCn1c(=O)n(C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=c1n(CF)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1CF', 'CCCn1c(=O)n(C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'C=Cn1c(=O)n(C=C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=C1N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)COCN1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'C=CN1CN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C1=O', 'COCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=C1NC(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)CN1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'O=c1n(CI)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'CCCCN1CN(CCCC)C(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C1=O', 'CCCCn1c(=O)n(CCC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCn1c(=O)[nH]c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCCn1c(=O)n(CCC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'O=C1NCN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'C=Cn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'C#CCn1c(=O)n(CC#C)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'N#CCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CCn1c(=O)n(CC)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'C#CCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'C=C1N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'O=c1n(CCBr)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'O=c1oc(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F', 'CC#CCn1c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c1=O', 'CN1CN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)CN(C)C1=O', 'CCCN1CN(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C(=O)N(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)C1=O', 'O=c1n(Cl)c(=O)n(C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F)c(=O)n1C(F)(F)C(F)(OC(F)(F)F)C(F)(F)F']

method = 'CAM-B3LYP'
basis = '6-31G(d,p)'
PCM = 'True'
EPS = 4.9
iteration = 3
debug = False

for i in range(len(name_list)):
    for j in range(iteration):
        utils.smile2xyz(name_list[i]+str(j)+'_0.xyz',smiles_list[i],randomSeed=None)

G0 = gaussian_cal.GaussianCal(method=method,basis=basis,charge='neu',wfn=True,debug=debug,PCM=PCM,EPS=EPS)
G1 = gaussian_cal.GaussianCal(method=method,basis=basis,charge='pos',wfn=True,debug=debug,PCM=PCM,EPS=EPS)

for i in range(len(name_list)):
    for j in range(iteration):
        if os.path.exists(f"{name_list[i]}{j}"):
            if debug:
                while utils.i8cpu_running():
                    time.sleep(60)
            G0.Run(f"{name_list[i]}{j}_0.xyz")
            if debug:
                while utils.i8cpu_running():
                    time.sleep(60)
            G1.Run(f"{name_list[i]}{j}_0.xyz")

