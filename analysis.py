'''
When you finsh the calculation.
You would like to get the result.
1. Every single EA from the files.

'''

import utils
import os

def IP_calculation(dir):
    csv = f'{dir}.csv'
    open(csv, 'w').close()
    try:
        if os.path.exists(f'{dir}/{dir}_2.log'):
            path1 = f'{dir}/{dir}_2.log'
            path2 = f'{dir}/{dir}_1.log'
        else:
            path1 = f'{dir}/{dir}_0.log'
            path2 = f'{dir}/{dir}_1.log'
        if utils.check_gaussian_log(path1) and utils.check_gaussian_log(path2):
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` > {1}".format(f'{path1}',csv))
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` >> {1}".format(f'{path2}',csv))
            with open(csv, 'r') as f:
                lines = f.readlines()
                cation = lines[-2].strip()
                neutral = lines[-1].strip()
                IP = (float(cation) - float(neutral)) * 27.2114
                os.system('rm {0}'.format(csv))
                os.system("echo  {0},{1} eV >> ip_val.csv".format(dir ,IP))
            # print (dir," cation energy(Ha):", cation, " neutral energy(Ha):", neutral, 'IP(eV):', IP)
                return(IP)
        else:
            os.system('rm {0}'.format(csv))
            return(0)
    except Exception as e:
        print(f"Error: {e}")
        os.system("echo  {0},Error eV >> ip_val.csv".format(dir))
        os.system('rm {0}'.format(csv))
        return(0)

def find_folders(folder):
    if os.path.exists(folder) and os.path.isdir(folder):
        subfolders = [f for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]
        return subfolders
    else:
        return []
    
    