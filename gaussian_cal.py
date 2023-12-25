import os
class GaussianCal:
    def __init__(self, method, basis, charge,opt = False,dispersion = False, PCM=False, wfn=True,debug=False):
        
        self.method = method
        self.basis = basis
        self.charge = charge

        self.opt = opt
        self.dispersion = dispersion
        self.PCM = PCM
        self.wfn = wfn
        self.debug = debug

    def Run(self,xyz_name):
        self.xyz_name = xyz_name
        dir_name = self.xyz_name.split('_')[0]
        # if not os.path.exists(dir_name):
        #     os.makedirs(dir_name)
        xyz_path = os.path.join(dir_name, self.xyz_name)
        header_path = os.path.join('header','gjf_header.txt')
        gjf_path = os.path.join(dir_name, self.xyz_name.replace(".xyz", ".gjf"))
        molename = self.xyz_name.split('.')[0]
        wfname = self.xyz_name.replace(".xyz", ".wfn")
        with open(header_path, 'r') as header_file:
            header_lines = header_file.read()

            modified_header_lines = header_lines.replace("__title__", molename)
            modified_header_lines = modified_header_lines.replace("__dir__", dir_name)
            # Parameters for Gaussian calculation
            modified_header_lines = modified_header_lines.replace("__method__", self.method)
            modified_header_lines = modified_header_lines.replace("__basis__", self.basis)

            charge_mapping = {'pos': "1 2", 'neg': "-1 2", 'neu': "0 1"}
            modified_header_lines = modified_header_lines.replace("__charge__", charge_mapping[self.charge])

            modified_header_lines = modified_header_lines.replace("__opt__", "opt" if self.opt else "")
            modified_header_lines = modified_header_lines.replace("__dispersion__", "em=gd3" if self.dispersion else "")
            modified_header_lines = modified_header_lines.replace("__PCM__", 'SCRF=(PCM,Solvent=Generic,Read)' if self.PCM else "")
            modified_header_lines = modified_header_lines.replace("__wfn__", f"out=wfn" if self.wfn else "")

        with open(xyz_path, 'r') as xyz_file:
            atomic_coordinates = xyz_file.readlines()[1:]  # skip the first line

        with open(gjf_path, 'w') as gjf_file:
            gjf_file.writelines(modified_header_lines)
            gjf_file.writelines(atomic_coordinates)
            gjf_file.write("\n")  # add a blank line
            if self.PCM:
                gjf_file.write("eps=2.05\n")
                gjf_file.write("\n")
            if self.wfn:
                gjf_file.write(f"{dir_name}/{wfname}\n")
            gjf_file.write("\n" * 4)  # add 4 blank lines, this is important for gaussian calculation

        sh_path = gjf_path.replace(".gjf", ".sh")
        sh_header_path = os.path.join('header', "sub_g09.txt")
        with open(sh_header_path, 'r') as file:
            file_content = file.read()
            modified_content = file_content.replace("__dir__", dir_name)
            modified_content = modified_content.replace("__name__", gjf_path.replace(".gjf", ""))
            modified_content = modified_content.replace("__changegjf__", gjf_path)
            modified_content = modified_content.replace("__changelog__", gjf_path.replace(".gjf", ".log"))
            modified_content = modified_content.replace("__debug__", 'i8cpu' if self.debug else 'F1cpu')
        with open(sh_path, 'w') as file:
            file.write(modified_content)
        os.system("sbatch {}".format(sh_path))
