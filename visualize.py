
import pandas as pd

df = pd.read_csv('Pubchem_FFKM.csv')
smiles = list(df['SMILES'])
cids = list(df['Name'])
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from tqdm import tqdm
def multiplot(smiles_list, cids,name):
    # 创建一个10x10的子图布局
    num_rows = len(smiles_list) // 10 + 1
    num_cols = 10
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols*5, num_rows*5))

    # 遍历SMILES表达式并在子图中显示
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            ax = axs[i // num_cols, i % num_cols]
            img = Draw.MolToImage(mol, size=(450, 450))
            ax.imshow(img)
            ax.axis('off')  # 关闭坐标轴
            ax.set_title(f"{cids[i]}", fontsize=40)  # 添加子图标题
    # for i in range(20, num_rows * num_cols):
    for i in range(len(smiles_list), num_rows * num_cols):
        axs.flatten()[i].axis('off')
    # 调整子图布局
    plt.tight_layout()
    plt.savefig(name,dpi=300,format='pdf')
    # 显示图像
    # plt.show()
for i in tqdm(range(0, len(smiles)//100 +1)):
    multiplot(smiles[i*100:(i+1)*100], cids[i*100:(i+1)*100],f'FFKM_{i}.pdf')
