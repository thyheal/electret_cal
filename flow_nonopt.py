import utils
import gaussian_cal
import pandas as pd

name_list = ['S','A','EHOPA','AEPYridine','AEP','APN','DBE','DIPEDA','OA','mXD']
smiles_list = [
'FC1(C(F)(C(F)(F)F)F)C(OC(F)(F)C(F)1F)(C(F)(C(F)(F)F)F)F',
'FC1(C(F)(C(F)(F)F)F)C(OC(F)(F)C(F)1F)(C(F)(C(O)=O)F)F',
'CCCCC(CC)COCCCNC(C(F)(F)C1(F)OC(F)(F)C(F)(F)C1(C(F)(C(O)=O)F)F)=O',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCC2=CN=CC=C2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN2CCCCC2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NC2=CC(C#N)=CC(C#N)=C2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN(CCCC)CCCC)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN(C(C)C)C(C)C)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCCCCCCC)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCC2=CC(CN)=CC=C2)=O)F)O1',
]


G0 = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-311G(d,p)',charge='pos',wfn=True,debug=True)
G0.log()

for i in range(len(name_list)):
    utils.smile2xyz(f'{name_list[i]}_0.xyz',smiles_list[i])
    utils.xyzcheck(f'{name_list[i]}_0.xyz',smiles_list[i])
    G0.Run(f'{name_list[i]}_0.xyz')
        # os.system(f'rm -r {name_list[i]}')


G1 = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-311G(d,p)',charge='neu',wfn=True,debug=True)
for i in range(len(name_list)):
    #提取上次高斯之后的坐标
    #检查经过一次高斯计算后是否还是原来的分子
    utils.xyzcheck(f'{name_list[i]}_0.xyz',smiles_list[i])
    G1.Run(f'{name_list[i]}_0.xyz')

G2 = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-311G(d,p)',charge='neg',wfn=True,debug=True)
for i in range(len(name_list)):
    #提取上次高斯之后的坐标
    #检查经过一次高斯计算后是否还是原来的分子
    utils.xyzcheck(f'{name_list[i]}_0.xyz',smiles_list[i])
    G2.Run(f'{name_list[i]}_0.xyz')

IP = []
IP_dict = {}
for _ in name_list:
    IP.append(float("{:.3f}".format((utils.IP_calculation(_)))))
table = pd.DataFrame({'Name': name_list, 'IP': IP})
table

EA = []
EA_dict = {}
for _ in name_list:
    EA.append(float("{:.3f}".format((utils.EA_calculation(_)))))
table = pd.DataFrame({'Name': name_list, 'EA': EA})
table
