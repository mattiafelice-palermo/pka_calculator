import os
import shutil
from rdkit import Chem

verbose = 0
dry_run = 0


class pKaCalculator:
    def __init__(self, molecule, charge, spin):

        self.molecule = molecule
        self.charge = charge
        self.spin = spin

    def single_point(self, xyz_dir, foldername, charge, spin):
        if not os.path.exists(foldername):
            os.makedirs(foldername)
        shutil.copyfile(
            xyz_dir + "/" + self.molecule + ".xyz", foldername + self.molecule + ".xyz"
        )

        os.chdir(foldername)

        if dry_run != 1:
            # os.system("echo $OMP_NUM_THREADS")
            os.system(
                f"xtb {self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water  > {self.molecule}.out 2>> {self.molecule}.out"
            )
        mol_file = [f for f in os.listdir(".") if f.endswith(".mol")][-1]

        start_mol = Chem.MolFromMolFile(
            mol_file, sanitize=False, removeHs=False, strictParsing=False
        )
        start_smiles = Chem.MolToSmiles(start_mol)

        os.chdir("../../")

        return start_smiles

    def optimization(self, xyz_dir, foldername, charge, spin):

        if not os.path.exists(foldername):
            os.makedirs(foldername)
        shutil.copyfile(
            xyz_dir + "/" + self.molecule + ".xyz", foldername + self.molecule + ".xyz"
        )

        os.chdir(foldername)

        if dry_run != 1:
            os.system(
                f"xtb {self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --ohess  > {self.molecule}.out 2>> {self.molecule}.out"
            )
        mol_file = [f for f in os.listdir(".") if f.endswith(".mol")][-1]

        end_mol = Chem.MolFromMolFile(
            mol_file, sanitize=False, removeHs=False, strictParsing=False
        )
        end_smiles = Chem.MolToSmiles(end_mol)

        for line in open(self.molecule + ".out", "r"):
            if "TOTAL FREE ENERGY" in line:
                energy = float(line.split()[-3])
        os.chdir("../../")

        return end_smiles, energy

    def deprotonate(self, xyz_dir, foldername, charge, spin):

        foldername = "deprotonate" + "/" + self.molecule + "/"

        if not os.path.exists("deprotonate/xyz_files"):
            os.makedirs("deprotonate/xyz_files")
        if not os.path.exists(foldername):
            os.makedirs(foldername)
        shutil.copyfile(
            xyz_dir + "/" + self.molecule + ".xyz", foldername + self.molecule + ".xyz"
        )

        os.chdir(foldername)

        if dry_run != 1:
            os.system(
                f"xtb {self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --ohess  > {self.molecule}_opt.out  2>> {self.molecule}_opt.out"
            )
            os.system(
                f"crest xtbopt.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -deprotonate > {self.molecule}_deprot.out  2>> {self.molecule}_deprot.out"
            )
        shutil.copyfile("deprotonated.xyz", "../xyz_files/" + self.molecule + ".xyz")

        os.chdir("../../")

    def compare_smiles(self, start_smiles, end_smiles, mol_type):

        with open("warnings.out", "a") as out:

            if verbose == 1:
                out.write(
                    f"\nMolecule: {self.molecule}\n {mol_type} Start SMILES: {start_smiles}\nEnd SMILES:   {end_smiles}\n"
                )
                if start_smiles != end_smiles:
                    out.write(
                        f"WARNING! Topology for molecule {self.molecule} {mol_type} has changed!\n"
                    )
            if "." in end_smiles:
                out.write(
                    f"WARNING!!! Molecule {self.molecule} {mol_type} has undergone dissociation!!!\n"
                )

    def calculate_pka(self):

        # radical cation starting single point
        start_smiles = self.single_point(
            "./xyz_files",
            "xtb_sp/" + self.molecule + "/",
            self.charge + 1,
            self.spin + 1,
        )
        # radical cation optimization
        end_smiles, protonated_energy = self.optimization(
            "./xyz_files",
            "xtb_opt/" + self.molecule + "/",
            self.charge + 1,
            self.spin + 1,
        )
        self.compare_smiles(start_smiles, end_smiles, "radical cation")

        # optimizing and deprotonating neutral molecule
        self.deprotonate(
            "./xyz_files", "deprotonate/" + self.molecule + "/", self.charge, self.spin
        )
        # neutral radical single point
        start_smiles = self.single_point(
            "./deprotonate/xyz_files",
            "xtb_sp_deprot/" + self.molecule + "/",
            self.charge,
            self.spin + 1,
        )
        # neutral radical optimization
        end_smiles, deprotonated_energy = self.optimization(
            "./deprotonate/xyz_files",
            "xtb_opt_deprot/" + self.molecule + "/",
            self.charge,
            self.spin + 1,
        )
        self.compare_smiles(start_smiles, end_smiles, "neutral radical")

        pka_correction = 164.22  # kcal/mol

        pka = -(
            (protonated_energy - deprotonated_energy) * 627.5 + 270.29 - pka_correction
        ) / (2.303 * 1.98720425864083 / 1000 * 298.15)

        # print(f'pKa for molecule {self.molecule} radical cation = {pka}')

        return self.molecule, pka
