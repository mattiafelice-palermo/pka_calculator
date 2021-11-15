import os, shutil
from multiprocessing import Process, Queue
from rdkit import Chem

xyz_dir = './xyz_files'
max_cores = 4
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

        self.molecule = molecule
        self.charge = charge
        self.spin = spin
        
        foldername = 'xtb_sp' + '/' + molecule + '/'

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + molecule + '.xyz', 
            foldername + molecule + '.xyz'
        )
        
        os.chdir(foldername)
        os.system(f'xtb {molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -P {cores_per_process} > {molecule}.out')
        
        mol_file = [ f for f in os.listdir('.') if f.endswith('.mol') ][-1]

        mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
        start_smiles = Chem.MolToSmiles(mol)

        print(f'{molecule} starting geometry: {start_smiles}')

        cores_list.get()
        
        return start_smiles


    def optimization(self):

        self.molecule = molecule
        self.charge = charge
        self.spin = spin

        foldername = 'xtb_opt' + '/' + molecule + '/'

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + molecule +'.xyz', 
            foldername + molecule +'.xyz'
        )
        
        os.chdir(foldername)
        os.system(f'xtb {molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --ohess -P {cores_per_process} > {molecule}.out')

        mol_file = [ f for f in os.listdir('.') if f.endswith('.mol') ][-1]

        mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
        end_smiles = Chem.MolToSmiles(mol)

        print(f'{molecule} final geometry: {end_smiles}')

        cores_list.get()
        
        return end_smiles


    def deprotonate(self, molecule, charge, spin):
        self.molecule = molecule
        self.charge = charge
        self.spin = spin

        foldername = 'deprotonate' + '/' + molecule + '/'

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + molecule +'.xyz', 
            foldername + molecule +'.xyz'
        )
        
        os.chdir(foldername)
        os.system(f'crest {molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -deprotonate -T {cores_per_process} > {molecule}.out')

        cores_list.get()    


    def compare_smiles(start_smiles, end_smiles):
        print(f'Start SMILES: {start_smiles}\n End SMILES: {end_smiles}')

    start_smiles = single_point()
    end_smiles = optimization()
    compare_smiles(start_smiles, end_smiles)

while molecule_list.empty() is False:
    if used_cores()+cores_per_process <= max_cores:
        try:
            molecule = molecule_list.get_nowait()
            charge = 0
            spin = 0
            worker = Process(target=Calculate_pka, args=(molecule, charge, spin))
            processes.append(worker)
            worker.start()
            cores_list.put(cores_per_process)
        except:
            print('Empty queue')

for process in processes:
    process.join()
