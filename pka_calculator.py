import os, shutil
from multiprocessing import Process, Queue
from rdkit import Chem

verbose = 1
dry_run = 1

max_cores = 4
cores_per_process = 1

molecule_list = Queue()
cores_list = Queue()

for xyz in os.listdir('./xyz_files'):
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


    def single_point(self, xyz_dir, foldername, charge, spin):
        
        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + self.molecule + '.xyz', 
            foldername + self.molecule + '.xyz'
        )
        
        os.chdir(foldername)

        if dry_run != 1:
            os.system(f'xtb {self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -P {cores_per_process} > {self.molecule}.out')
        
        mol_file = [ f for f in os.listdir('.') if f.endswith('.mol') ][-1]

        start_mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
        start_smiles = Chem.MolToSmiles(start_mol)

        os.chdir('../../')
        
        return start_smiles


    def optimization(self, xyz_dir, foldername, charge, spin):

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + self.molecule +'.xyz', 
            foldername + self.molecule +'.xyz'
        )
        
        os.chdir(foldername)
    
        if dry_run != 1:
            os.system(f'xtb {self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --ohess -P {cores_per_process} > {self.molecule}.out')

        mol_file = [ f for f in os.listdir('.') if f.endswith('.mol') ][-1]

        end_mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
        end_smiles = Chem.MolToSmiles(end_mol)

        for line in open(self.molecule+'.out', 'r'):
            if "TOTAL FREE ENERGY" in line:
                energy = float(line.split()[-3])
        
        os.chdir('../../')
        
        return end_smiles, energy


    def deprotonate(self, xyz_dir, foldername, charge, spin):

        foldername = 'deprotonate' + '/' + self.molecule + '/'

        if not os.path.exists('deprotonate/xyz_files'):
            os.makedirs('deprotonate/xyz_files')

        if not os.path.exists(foldername):
            os.makedirs(foldername)
            
        shutil.copyfile(
            xyz_dir + '/' + self.molecule +'.xyz', 
            foldername + self.molecule +'.xyz'
        )
        
        os.chdir(foldername)

        if dry_run != 1:
            os.system(f'crest {self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -deprotonate -T {cores_per_process} > {molecule}.out')
        
        shutil.copyfile(
            'deprotonated.xyz', 
            '../xyz_files/' + self.molecule + '.xyz'
        )
        
        os.chdir('../../')

    def compare_smiles(self, start_smiles, end_smiles):

        with open('logfile.out', 'a') as out:

            if verbose == 1:
                out.write(f'\nMolecule: {self.molecule}\nStart SMILES: {start_smiles}\nEnd SMILES:   {end_smiles}\n')
                if start_smiles != end_smiles:
                    out.write(f'WARNING! Topology for molecule {self.molecule} has changed!\n')

            if '.' in end_smiles:
                out.write(f'WARNING!!! Molecule {self.molecule} has undergone dissociation!!!\n')


    def calculate_pka(self):

        start_smiles = self.single_point(
            './xyz_files', 
            'xtb_sp/'+self.molecule+'/', 
            self.charge, 
            self.spin
        )
        end_smiles, protonated_energy = self.optimization(
            './xyz_files', 
            'xtb_opt/'+self.molecule+'/', 
            self.charge, 
            self.spin
        )
        self.compare_smiles(start_smiles, end_smiles)

        self.deprotonate(
            './xyz_files', 
            'deprotonate/'+self.molecule+'/',
            self.charge, 
            self.spin
        )

        start_smiles = self.single_point(
            './deprotonate/xyz_files', 
            'xtb_sp_deprot/'+self.molecule+'/', 
            self.charge-1, 
            self.spin
        )
        end_smiles, deprotonated_energy = self.optimization(
            './deprotonate/xyz_files', 
            'xtb_opt_deprot/'+self.molecule+'/', 
            self.charge-1, 
            self.spin
        )
        self.compare_smiles(start_smiles, end_smiles)

        pka_correction = 164.22  # kcal/mol

        pka = ((protonated_energy - deprotonated_energy) * 627.5 + 270.29 - pka_correction) \
            / (2.303 * 1.98720425864083 / 1000 * 298.15) 

        print(f'pKa for molecule {self.molecule} = {pka}')

        cores_list.get() 

        return pka

with open('logfile.out', 'w') as out:
    out.write('--- Running pKa calculation ---\n')

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

