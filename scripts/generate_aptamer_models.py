#! /usr/bin/env python
import os
import sys
import time
import argparse
import pandas as pd
sys.path.insert(1, './scripts')
import MutateRna
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager

###############
# Parses script accepts RNA and DNA aptamers sequences with folds in dot-bracket notation and produces pdb models for docking
###############

parser = argparse.ArgumentParser(description='Prepare aptamer 3D pdb models and ligands for docking')
parser.add_argument('-i', nargs=1, help='Input csv file with aptamer sequences, folding data and ligands', type=str, required=True)
parser.add_argument('-f', '--fold', nargs=1, help='What secondary structure prediction to use. Options: seqfold, rnafold', type=str, required=True)
args = parser.parse_args()

if not args.fold[0] in ['seqfold', 'rnafold']:
   raise Exception('A wrong option provided with --fold flag. Use one of these: seqfold, rnafold')

### read and inspect records from csv file
print(f"Reading {args.i[0]} file...\n")
aptamer_data = pd.read_csv(args.i[0])
rna_count, dna_count, discarded = 0, 0, 0
aptamers = []
renamed_aptamers = dict()
# inspect if the parameters are correct
for i, (name, seq, fold) in enumerate(zip(aptamer_data['aptamer_name'], aptamer_data['sequence'], aptamer_data[args.fold[0]])):
    # check the name length. RNAComposer requires name no longer than 15 characters
    if len(name) > 15:
        renamed_aptamers[name[-15:]] = name
        name = name[-15:]
        
    # check that the sequence is either DNA or RNA
    if 'U' in seq:
        if set(seq) > set('ACGU'):
            print(f'Skipping sequence at the line {i+1}. The sequence is neither DNA or RNA: {seq}')
            discarded += 1
            continue
        else:
            is_dna = False
            rna_count += 1
    else:
        if set(seq) > set('ACGT'):
            print(f'Skipping sequence at the line {i+1}. The sequence is neither DNA or RNA: {seq}')
            discarded += 1
            continue
        else:
            is_dna = True
            dna_count += 1
            
    if len(set(fold) - set('.()[]{}')) > 0:
        print(f'Skipping sequence at the line {i+1}. The dot-bracket notation contains unknown symbols: {fold}')
        discarded += 1
        continue
    if len(fold) != len(seq):
        print(f'Skipping the record {name}, the sequence and the fold are of different length: {seq}, {fold}')
        discarded += 1
        continue
        
    aptamers.append([name, seq, fold, is_dna])

del aptamer_data

print("Finished loading the aptamers")
print(f"Collected {len(aptamers)} records")
print(f"{dna_count} DNA aptamers")
print(f"{rna_count} RNA aptamers")
print(f"Dropped {discarded} records from original file\n")

### Request pdb structures from RNAComposer server

print("Requesting pdb structures from RNAComposer server")
os.makedirs(f"2_docking/aptamers_pdb", exist_ok=True)

#set a directory for download
options = webdriver.ChromeOptions()
options.add_argument('headless')
prefs = {"download.default_directory": f"{os.getcwd()}/2_docking/aptamers_pdb"}
options.add_experimental_option("prefs", prefs);

# create a driver to work with
driver = webdriver.Chrome(ChromeDriverManager().install(), options = options)

print("\nBeginning to request the aptamer models...")
for aptamer in aptamers:
    name, seq, fold, is_dna = aptamer
    
    print(f"Requesting the model for {name}")
    # Open the RNA_Composer web-server
    driver.get('https://rnacomposer.cs.put.poznan.pl')
    # Select the input box by id
    input_box = driver.find_element_by_id('input')
    # clear the input field
    input_box.clear()
    # Send input
    input_box.send_keys(f'>{name}\n')
    input_box.send_keys(f'{seq.translate(str.maketrans("T", "U"))}\n')
    input_box.send_keys(fold)

    # Find 'compose' button
    send_button = driver.find_element_by_name('send')
    # Click to send the request
    send_button.click()

    # wait until the result is loaded
    result = None
    while result == None:
        try:
            result = driver.find_element_by_link_text(f'{name}.pdb')
        except:
            time.sleep(1)  # wait more for a result to be ready

    result.click()
    time.sleep(5)  # wait for a file to be downloaded

    print(f"Downloaded RNAComposer model for an aptamer {name}")

### Convert DNA aptamers from RNA if needed
print(f"Converting {dna_count} pdb structures from RNA to DNA")
for aptamer in aptamers:
    name, seq, fold, is_dna = aptamer
    if is_dna:
        MutateRna.convert_pdb(f"{os.getcwd()}/2_docking/aptamers_pdb/{name}.pdb")
        print(f"Converted {name}.pdb from RNA to DNA")
        
    # rename the model if needed
    if name in renamed_aptamers.keys():
        os.rename(f'{os.getcwd()}/2_docking/aptamers_pdb/{name}.pdb',
                  f'{os.getcwd()}/2_docking/aptamers_pdb/{renamed_aptamers[name]}.pdb')
