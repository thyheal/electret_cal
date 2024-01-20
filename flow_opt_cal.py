import utils
import gaussian_cal
import pandas as pd

smiles_list = ['C','N','CC','NN']
name_list = ['Ch4','NH3','Ethane','Hydrazine']

G0 = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31+G(d,p)', charge='pos',wfn=True,debug=True, opt=True)
for i in range(len(name_list)):
    utils.smile2xyz(f'{name_list[i]}_0.xyz',smiles_list[i])
    G0.Run(f'{name_list[i]}_0.xyz')

G1 = gaussian_cal.GaussianCal(method='cam-B3LYP',basis='6-31+G(d,p)',charge='neu',wfn=True,debug=True)
for i in range(len(name_list)):
    #提取上次高斯之后的坐标
    utils.log2xyz(f'{name_list[i]}_0.log')
    #检查经过一次高斯计算后是否还是原来的分子
    if utils.xyzcheck(f'{name_list[i]}_1.xyz',smiles_list[i]):
        G1.Run(f'{name_list[i]}_1.xyz')
    else:
        print(name_list[i],'_1.log is wrong.')

IP = []
IP_dict = {}
for _ in name_list:
    IP.append(float("{:.3f}".format((utils.IP_calculation(_)))))
table = pd.DataFrame({'Name': name_list, 'IP': IP})

        
