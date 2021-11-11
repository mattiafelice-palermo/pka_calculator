import os, shutil
from multiprocessing import Process, Queue

xyz_dir = './xyz_files'

molecule_list = Queue()
core_list = Queue()

max_cores = 8
cores_per_process = 4

for xyz in os.listdir(xyz_dir):
    molecule_list.put(xyz.strip('.xyz'))
    print(xyz)

print(molecule_list)

def used_cores():
    ncores = cores_per_process*core_list.qsize()
    return ncores


def B97_opt(molecule, charge, spin):
    
    foldername = 'B97_opt' + '/' + molecule + '/'

    if not os.path.exists(foldername):
        os.makedirs(foldername)

    shutil.copyfile(
        xyz_dir + '/' + molecule +'.xyz', 
        foldername + molecule +'.xyz'
    )

    with open(foldername + molecule + '.inp', 'w') as file:

        file.write(f"""%pal
  nprocs 1
end
%maxcore 375

! B97-3c Opt Freq

xyzfile {charge} {spin} {molecule}.xyz
""")


def xtb_opt(molecule, charge, spin):
    
    foldername = 'xtb_opt' + '/' + molecule + '/'

    if not os.path.exists(foldername):
        os.makedirs(foldername)
        
    shutil.copyfile(
        xyz_dir + '/' + molecule +'.xyz', 
        foldername + molecule +'.xyz'
    )

    os.system(f'xtb {molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -P {cores_per_process} > {molecule}.out')
    core_list.get()

calculation = xtb_opt
charge = 0
spin = 0

processes = []

while molecule_list.empty() is False:
    if used_cores()+cores_per_process <= max_cores:
        try:
            molecule = molecule_list.get_nowait()
            worker = Process(target=xtb_opt, args=(molecule, charge, spin))
            processes.append(worker)
            worker.start()
            core_list.put(cores_per_process)
        except:
            print('Empty queue')

for process in processes:
    process.join()
    