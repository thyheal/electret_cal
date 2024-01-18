import utils
import gaussian_cal
import pandas as pd

name_list = ['Benzoqionone',
             'TCNQ',
             'F4TCNQ']
smiles_list = ['C1=CC(=O)C=CC1=O',
               'C1=CC(=C(C#N)C#N)C=CC1=C(C#N)C#N',
               'C(#N)C(=C1C(=C(C(=C(C#N)C#N)C(=C1F)F)F)F)C#N']

G0 = gaussian_cal.GaussianCal(method='B3LYP',basis='6-31G',charge='pos',wfn=True,debug=True)
for i in range(len(name_list)):
    if utils.xyzcheck(f'{name_list[i]}_0.xyz',smiles_list[i]):
        G0.Run(f'{name_list[i]}_0.xyz')
    else:
        print(name_list[i],"_0.xyz is wrong.")
        # os.system(f'rm -r {name_list[i]}')


G1 = gaussian_cal.GaussianCal(method='B3LYP',basis='6-31G',charge='neu',wfn=True,debug=True)
for i in range(len(name_list)):
    #提取上次高斯之后的坐标
    utils.log2xyz(f'{name_list[i]}_0.log')
    #检查经过一次高斯计算后是否还是原来的分子
    if utils.xyzcheck(f'{name_list[i]}_1.xyz',smiles_list[i]):
        G1.Run(f'{name_list[i]}_1.xyz')
    else:
        print(name_list[i],"_1.log is wrong.")

IP = []
IP_dict = {}
for _ in name_list:
    IP.append(float("{:.3f}".format((utils.IP_calculation(_)))))
table = pd.DataFrame({'Name': name_list, 'IP': IP})
table
        
