import matplotlib.pyplot as plt
from rdkit.Chem import Draw
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Draw
from rdkit import Chem
from datetime import datetime
import matplotlib.pyplot as plt

def multi_mol_plot(smiles_list, name_list, prop=False, save_name=False, cols=10):
    current_time = datetime.now()
    time_string = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    file_name = f"output_{time_string}.txt"
    if not save_name:
        save_name = file_name
        
    # 创建一个10x10的子图布局
    num_rows = (len(smiles_list) + cols - 1) // cols
    num_cols = cols
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols*5, num_rows*5))
    # 遍历SMILES表达式并在子图中显示
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            ax = axs[i // num_cols, i % num_cols]
            img = Draw.MolToImage(mol, size=(450, 450))
            ax.imshow(img)
            ax.axis('off')  # 关闭坐标轴
            if prop:
                ax.set_title(f"{name_list[i], prop[i]}", fontsize=30)
            else:
                ax.set_title(f"{name_list[i]}", fontsize=30)

    # 隐藏多余的子图
    for i in range(len(smiles_list), num_rows * num_cols):
        axs.flatten()[i].axis('off')

    # 调整子图布局
    plt.tight_layout()

    # 保存图像
    plt.savefig(save_name, dpi=300, format='pdf')

