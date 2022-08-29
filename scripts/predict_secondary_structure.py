#! /usr/bin/env python
import pandas as pd
import argparse
import os
from tqdm import tqdm
import subprocess
from subprocess import Popen, PIPE, run, DEVNULL, STDOUT

###############
# This script automates secondary structure prediction for aptamer sequences
###############

parser = argparse.ArgumentParser(description='Run secondary structure prediction for aptamers')
# input as fasta
parser.add_argument('-f', '--fasta', nargs=1, help='Fasta file with input sequences', type=str, required=False)
parser.add_argument('-l', '--ligand_name', nargs=1, help='Ligand molecule name', type=str, required=False)
parser.add_argument('-s', '--ligand_smiles', nargs=1, help='Ligand molecule as SMILES string (e.x. "CCC")', type=str, required=False)

# input as csv table
parser.add_argument('-c', '--csv', nargs=1, help='CSV table with aptamer sequences and ligands', type=str, required=False)

parser.add_argument('--seqfold', help='Use seqfold', action='store_true', required=False)
parser.add_argument('--rnafold', help='Use rnafold', action='store_true', required=False)
parser.add_argument('--nupack', help='Use nupack', action='store_true', required=False)


args = parser.parse_args()

os.makedirs('1_secondary_structure', exist_ok=True)

### check and parse the input
if not (args.seqfold or args.rnafold or args.nupack):
    raise Exception(f"No secondary structure prediction method specified. Provide at least one")

if args.fasta is None:
    if args.csv is None:
        raise Exception("No input provided. Provide either fasta file or csv table")
    else:
        aptamers = pd.read_csv(args.csv[0])[['aptamer_name', 'sequence', 'ligand_name', 'ligand_smiles']]
        def crop(seq):
            if len(seq) > 15:
                seq = seq[:15]
            return seq
        aptamers['aptamer_name'] = aptamers.apply(lambda row: crop(row.aptamer_name), axis=1)
        # alse prepare a fasta file
        with open('1_secondary_structure/aptamers.fasta', 'w') as f:
            for name, seq in zip(aptamers['aptamer_name'], aptamers['sequence']):
                f.write(f">{name}\n")
                f.write(f"{seq}\n")
        fasta_file = '1_secondary_structure/aptamers.fasta'
        print(f"Loaded data on {aptamers.shape[0]} with their ligands from csv\n")
else:
    if args.ligand_name is None or args.ligand_smiles is None:
        raise Exception("Provide both ligand name and ligand smiles")
    # parse fasta
    fasta_file = args.fasta[0]
    with open(args.fasta[0], 'r') as f:
        lines = [line.strip() for line in f]
    aptamers = []
    for i in range(0, len(lines), 2):
        name, seq = lines[i], lines[i+1]
        if not name.startswith('>'):
            raise Exception(f"Error parsing fasta file on line {i+1}: name. The line does not appear to be a correct fasta name")
        else:
            name = name[1:]
            if len(name)>15:
                name = name[-15:]
            aptamers.append([name, seq, args.ligand_name[0], args.ligand_smiles[0]])
            
    aptamers = pd.DataFrame(data=aptamers, columns=['aptamer_name', 'sequence', 'ligand_name', 'ligand_smiles'])
    print(f"Loaded data on {aptamers.shape[0]} aptamers from fasta\n")
    
### Run secondary structure prediction

### seqfold
if args.seqfold:

    folds, dGs = [], []
    for seq in tqdm(aptamers['sequence'], desc="Running seqfold secondary structure prediction:"):
        seqfold_args = ['seqfold', '-d', seq]
        result = run(seqfold_args, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        output = result.stdout.split("\n")
        folds.append(output[1])
        dGs.append(float(output[2]))
    
    aptamers['seqfold'] = folds
    aptamers['seqfold_ddg'] = dGs
    
    print(f"Finished secondary structure prediction with seqfold!\n")
    
    
### rnafold
if args.rnafold:

    print("Running rnafold secondary structure prediction. The results will be written to 1_secondary_structure/rnafold_results.txt")
    rnafold_args = ['rnafold', '-p', '--MEA', '-g', f"{os.getcwd()}/{fasta_file}"]
    os.makedirs('1_secondary_structure/rnafold_outputs', exist_ok=True)
    os.chdir('1_secondary_structure/rnafold_outputs')
    out_file = open('../rnafold_results.txt', 'w')
    subprocess.call(rnafold_args, stdout=out_file)
    os.chdir(f"../../")
    out_file.close()
    
    with open('1_secondary_structure/rnafold_results.txt', 'r') as f:
        lines = [line.strip() for line in f]
    folds, dGs = [], []
    for i in range(0, len(lines), 7):
        dat = lines[i+2].split()
        folds.append(dat[0])
        dGs.append(float(dat[-1].translate({ord('('):None, ord(')'): None})))

    aptamers['rnafold'] = folds
    aptamers['rnafold_ddg'] = dGs
    print(f"Finished secondary structure prediction with rnafold!\n")
    
### nupack
if args.nupack:
    from nupack import *
    
    print("Running nupack secondary structure prediction...")
    
    folds, dGs = [], []
    model_dna = Model(material='dna', celsius=25, sodium=0.120, magnesium=0.005)
    for seq in tqdm(aptamers['sequence'], desc="Running nupack secondary structure prediction:"):

        A = Strand(seq, name='A')
        
        # specify complex set, we only study folding
        set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1))
        
        # calculate the partition function for each complex in the complex set
        complex_results1 = complex_analysis(complexes=set1, model=model_dna, compute=['pfunc', 'pairs', 'mfe'])
        
        # return MFE structure and energy
        for key in complex_results1.complexes.keys():
            fold, ddg = complex_results1.complexes[key].mfe[0].structure, float(complex_results1.complexes[key].free_energy)

        folds.append(fold)
        dGs.append(ddg)
    
    aptamers['nupack'] = folds
    aptamers['nupack_ddg'] = dGs
    print(f"Finished secondary structure prediction with nupack!\n")
    
    
    
### save data on foldings
aptamers.to_csv('1_secondary_structure/folded_aptamers.csv')
    
        
        
        
    
