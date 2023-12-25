import utils
import gaussian_cal
import logging

name_list = ['A11','AEP5','APN33']
smiles_list = ['FC1(C(F)(C(F)(F)F)F)C(OC(F)(F)C(F)1F)(C(F)(C(O)=O)F)F',
            'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN2CCCCC2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NC2=CC(C#N)=CC(C#N)=C2)=O)F)O1',]

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
        
