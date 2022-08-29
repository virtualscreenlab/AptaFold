#! /usr/bin/env python
import os
import glob
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from subprocess import Popen, PIPE, run, DEVNULL, STDOUT

###############
# This script automates small molecule docking to the aptamers
###############

parser = argparse.ArgumentParser(description='Run aptamer-target docking')
parser.add_argument('-a', '--aptamers', nargs=1, help='Directory containing aptamer pdb files for docking', type=str, required=True)
parser.add_argument('-f', '--flags', nargs=1, help='flags file for docking', type=str, required=True)
parser.add_argument('-l', '--ligand_data', nargs=1, help='A csv file with aptamer names and ligand data', type=str, required=True)
args = parser.parse_args()

### list aptamers for docking and convert them to pdbqt
print("Started aptamer conversion to pdbqt")
aptamers = glob.glob(f"{args.aptamers[0]}/*.pdb")
os.makedirs('aptamers_pdbqt', exist_ok=True)
for aptamer in aptamers:
    name = aptamer.split('/')[-1][:-4]
    process = Popen(["obabel", "-ipdb", aptamer, "-opdbqt", "-O", f"aptamers_pdbqt/{name}.pdbqt"], stdout=DEVNULL, stderr=STDOUT)
    returncode = process.wait()
    with open(f"aptamers_pdbqt/{name}.pdbqt", 'r') as f:
        lines = f.readlines()
    with open(f"aptamers_pdbqt/{name}.pdbqt", 'w') as f:
        for line in lines:
            if line.startswith('ATOM'):
                f.write(line)
    print(f"Converted {aptamer}")
print(f"\nConverted {len(aptamers)} aptamers to pdbqt files\n")

### list ligands to use for docking and convert them to pdbqt

aptamers_data = pd.read_csv(args.ligand_data[0])
ligands_list = defaultdict(list)
os.makedirs('ligands/pdbqt', exist_ok=True)
os.makedirs('ligands/mol2', exist_ok=True)
for aptamer_name, ligand_name, smiles in zip(aptamers_data['aptamer_name'], aptamers_data['ligand_name'], aptamers_data['ligand_smiles']):
    if not ligand_name in ligands_list.keys():
        process = Popen(["obabel", f'-:{smiles}', "-O", f"ligands/mol2/{ligand_name}.mol2", "--gen3d", "slow"], stdout=DEVNULL, stderr=STDOUT)
        returncode = process.wait()
        process = Popen(["obabel", "-imol2", f"ligands/mol2/{ligand_name}.mol2", "-opdbqt", "-O", f"ligands/pdbqt/{ligand_name}.pdbqt"], stdout=DEVNULL, stderr=STDOUT)
        returncode = process.wait()
        print(f"Converted {ligand_name} ligand to pdbqt")
    ligands_list[ligand_name].append(aptamer_name)
print(f"Found {len(ligands_list.keys())} ligands and converted them to pdbqt.\n")
    
### find aptamer pockets for docking

print(f"Starting search for binding pockets in aptamers\n")
all_binding_pockets = defaultdict(list)
for aptamer in aptamers:
    name = aptamer.split('/')[-1][:-4]
    # run fpocket
    process = Popen(["fpocket", "-f", aptamer], stdout=DEVNULL, stderr=STDOUT)
    returncode = process.wait()

    # parse pockets centroids
    pockets = defaultdict(list)
    with open(f"{args.aptamers[0]}/{name}_out/{name}_out.pdb", 'r') as f:
        for line in f:
            dat = line.split()
            if dat[2] == 'APOL':
                if len(dat) == 12:
                    x, y, z = float(dat[6]), float(dat[7]), float(dat[8])
                else:
                    if dat[6].count('.') > 1:
                        sep = dat[6][1:].find('-')
                        x, y, z = float(dat[6][:sep+1]), float(dat[6][sep+1:]), float(dat[7])
                    elif dat[7].count('.') > 1:
                        sep = dat[7][1:].find('-')
                        x, y, z = float(dat[6]), float(dat[7][:sep+1]), float(dat[7][sep+1:])
                pockets[dat[5]].append([x, y, z])
    print(f"Found {len(pockets)} distinct pockets in fpocket output for {name}.pdb")
    for pocket_id in pockets:
        pocket_atoms = np.array(pockets[pocket_id])
        all_binding_pockets[name].append(np.mean(pocket_atoms, axis = 0))
        
print(f"\nFinished search for binding pockets in aptamers\n")

### run docking

print(f"Starting docking runs...\n")

# read flags for docking
flags = []
with open(args.flags[0], 'r') as f:
    for line in f:
        flags.append(line.strip().split(' '))

# run docking for all ligands and their aptamers
for ligand in ligands_list:
    print(f"Running docking for {ligand} ligand")
    aptamers = ligands_list[ligand]
    os.makedirs('docking_out/pdbqt', exist_ok=True)
    os.makedirs('docking_out/pdb', exist_ok=True)
    docking_data = []
    for name in aptamers:
        print(f"Running docking for the aptamer {name}")

        for i, pocket in enumerate(all_binding_pockets[name]):
            print(f"Docking into the pocket {i}")
            docking_flags = [os.environ['VINA'],
                           "--receptor", f"aptamers_pdbqt/{name}.pdbqt",
                           "--ligand", f"ligands/pdbqt/{ligand}.pdbqt",
                           "--center_x", str(pocket[0]),
                           "--center_y", str(pocket[1]),
                           "--center_z", str(pocket[2])]
            for flag in flags:
                docking_flags.extend(flag)
            docking_flags.extend(["--out", f"docking_out/pdbqt/{name}_pocket{i}_{ligand}.pdbqt"])

            # run docking and parse output
            result = run(docking_flags, stdout=PIPE, stderr=PIPE, universal_newlines=True)
            output = result.stdout.split("\n")
            
            start_reading = False
            for line in output:
                if start_reading and not ("Writing output" in line) and (len(line) > 0):
                    dat = line.split()
                    docking_data.append([name, ligand, i, int(dat[0]), float(dat[1]), float(dat[2]), float(dat[3])])
                if "-----+------------+----------+----------" in line:
                    start_reading = True
            
            # convert pdbqt output to pdb
            process = Popen(["obabel", "-ipdbqt", f"docking_out/pdbqt/{name}_pocket{i}_{ligand}.pdbqt", "-opdb", "-O", f"docking_out/pdb/{name}_pocket{i}_{ligand}.pdb"], stdout=DEVNULL, stderr=STDOUT)
            returncode = process.wait()
            
        print(f"Finished docking for the aptamer {name}\n")
    
# save docking data
docking_data = pd.DataFrame(data=docking_data,
                            columns=['aptamer', 'ligand', 'pocket', 'pose', 'affinity (kcal/mol)', 'rmsd l.b.', 'rmsd u.b.'])
docking_data.to_csv('docking_data.csv')


