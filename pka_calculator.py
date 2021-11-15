import os, shutil
from multiprocessing import Process, Queue
from rdkit import Chem

xyz_dir = './xyz_files'
max_cores = 8
cores_per_process = 1

molecule_list = Queue()
cores_list = Queue()

for xyz in os.listdir(xyz_dir):
    molecule_list.put(xyz.strip('.xyz'))

processes = []

def used_cores():
    ncores = cores_per_process*cores_list.qsize()
    return ncores


class Calculate_pka:

    def __init__(self, molecule, charge, spin):

        self.molecule = molecule
        self.charge = charge
        self.spin = spin


    def single_point(self):
        
        foldername = 'xtb_sp' + '/' + self.molecule + '/'

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + self.molecule + '.xyz', 
            foldername + self.molecule + '.xyz'
        )
        
        os.chdir(foldername)
        os.system(f'xtb {self.molecule}.xyz --gfn2 --chrg {self.charge} --uhf {self.spin} --alpb water -P {cores_per_process} > {self.molecule}.out')
        
        mol_file = [ f for f in os.listdir('.') if f.endswith('.mol') ][-1]

        mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
        start_smiles = Chem.MolToSmiles(mol)

        os.chdir('../../')
        
        return start_smiles


    def optimization(self):

        foldername = 'xtb_opt' + '/' + self.molecule + '/'

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + self.molecule +'.xyz', 
            foldername + self.molecule +'.xyz'
        )
        
        os.chdir(foldername)
        os.system(f'xtb {self.molecule}.xyz --gfn2 --chrg {self.charge} --uhf {self.spin} --alpb water --ohess -P {cores_per_process} > {self.molecule}.out')

        mol_file = [ f for f in os.listdir('.') if f.endswith('.mol') ][-1]

        mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
        end_smiles = Chem.MolToSmiles(mol)

        os.chdir('../../')
        
        return end_smiles


    def deprotonate(self):

        foldername = 'deprotonate' + '/' + self.molecule + '/'

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + self.molecule +'.xyz', 
            foldername + self.molecule +'.xyz'
        )
        
        os.chdir(foldername)
        os.system(f'crest {self.molecule}.xyz --gfn2 --chrg {self.charge} --uhf {self.spin} --alpb water -deprotonate -T {cores_per_process} > {molecule}.out')
        os.chdir('../../')

    def compare_smiles(self, start, end):

        print(f'Molecule: {self.molecule}\nStart SMILES: {start}\nEnd SMILES:   {end}\n')
        
        if start != end:
            print(f'WARNING! Topology for molecule {self.molecule} has changed!\n')

    def calculate_pka(self):
        start = self.single_point()
        end = self.optimization()
        self.compare_smiles(start, end)
        self.deprotonate()
        cores_list.get() 



while molecule_list.empty() is False:
    if used_cores()+cores_per_process <= max_cores:
        try:
            molecule = molecule_list.get_nowait()
            singlet = Calculate_pka(molecule, 0, 0)
            worker = Process(target=singlet.calculate_pka, args=())
            processes.append(worker)
            worker.start()
            cores_list.put(cores_per_process)
        except:
            print('Empty queue')

for process in processes:
    process.join()

