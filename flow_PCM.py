# Used
import utils
import gaussian_cal
import pandas as pd

smiles_list = ['C1=CC(=O)C=CC1=O',
               'C1=CC(=C(C#N)C#N)C=CC1=C(C#N)C#N',
               'C(#N)C(=C1C(=C(C(=C(C#N)C#N)C(=C1F)F)F)F)C#N']
name_list = ['Benzoqionone',
             'TCNQ',
             'F4TCNQ']

# G0 = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31+G(d,p)',polar=True,Volume=True, charge='pos',wfn=True,debug=True)
# for i in range(len(name_list)):
#     # utils.smile2xyz(f'{name_list[i]}_0.xyz',smiles_list[i])
#     G0.Run(f'{name_list[i]}_0.xyz')

GB = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31G(d,p)',EPS =1.131, PCM=True,polar=True,Volume=True, opt=True, charge='pos',wfn=True,debug=True)
GT = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31G(d,p)',EPS =1.589, PCM=True,polar=True,Volume=True, opt=True,charge='pos',wfn=True,debug=True)
GF = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31G(d,p)',EPS =1.752, PCM=True,polar=True,Volume=True, opt=True,charge='pos',wfn=True,debug=True)
i = 0
GB.Run(f'{name_list[i]}_0.xyz')
i = 1
GT.Run(f'{name_list[i]}_0.xyz')
i = 2
GF.Run(f'{name_list[i]}_0.xyz')


GB = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31G(d,p)',EPS =1.131, PCM=True, charge='neu',wfn=True,debug=True)
GT = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31G(d,p)',EPS =1.589, PCM=True, charge='neu',wfn=True,debug=True)
GF = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31G(d,p)',EPS =1.752, PCM=True, charge='neu',wfn=True,debug=True)
for i in range(len(name_list)):
    #提取上次高斯之后的坐标
    utils.log2xyz(f'{name_list[i]}_0.log')
    #检查经过一次高斯计算后是否还是原来的分子
    utils.xyzcheck(f'{name_list[i]}_1.xyz',smiles_list[i])
i = 0
GB.Run(f'{name_list[i]}_1.xyz')
i = 1
GT.Run(f'{name_list[i]}_1.xyz')
i = 2
GF.Run(f'{name_list[i]}_1.xyz')


IP = []
for _ in name_list:
    IP.append(float("{:.3f}".format((utils.IP_calculation(_)))))

table = pd.DataFrame({'Name': name_list, 'IP': IP})
table.to_csv('IP.csv',index=False)
table

alpha = [45.345,392.626,471.264]
V = [4529.919,10015.856,9850.078]
for i in range(3):
    print(3/(1-4/3*alpha[i]/V[i]*3.142)-2)

