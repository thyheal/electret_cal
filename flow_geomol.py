import utils
import gaussian_cal
import pandas as pd

name_list = ['Benzoqionone',
             'TCNQ',
             ]
smiles_list = ['C1=CC(=O)C=CC1=O',
               'C1=CC(=C(C#N)C#N)C=CC1=C(C#N)C#N',
               ]

G0 = gaussian_cal.GaussianCal(method='B3LYP',basis='6-31G(d,p)',charge='pos',wfn=True,debug=True)
G0.log()

for i in range(len(name_list)):
    utils.smile2xyz(f'{name_list[i]}_0.xyz',smiles_list[i])
    utils.xyzcheck(f'{name_list[i]}_0.xyz',smiles_list[i])
    G0.Run(f'{name_list[i]}_0.xyz')
        # os.system(f'rm -r {name_list[i]}')


G1 = gaussian_cal.GaussianCal(method='B3LYP',basis='6-31G(d,p)',charge='neg',wfn=True,debug=True)
for i in range(len(name_list)):
    #提取上次高斯之后的坐标
    #检查经过一次高斯计算后是否还是原来的分子
    utils.xyzcheck(f'{name_list[i]}_0.xyz',smiles_list[i])
    G1.Run(f'{name_list[i]}_0.xyz')



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
