# AptaFold workflow
AptaFold is a python3 based workflow for small organic molecule docking with DNA and RNA aptamers starting from oligonucleotide sequence. 

## Setup

Before cloning this repository and starting the workflow, setup the environment and dependecies:
1. It is highly recommended to create a separate python environment for AptaFold jobs. `Conda` can be convenient in this case as it provides easy access to some necessary tools and packages:
```
conda create --name myenv
conda activate myenv
```
2. Install requred python libraries with conda or pip. **Note:** when using pip, make sure you are using it from the current conda environment by running `which pip`.
```
conda install -c anaconda pandas numpy
conda install -c conda-forge selenium biopython tqdm

#or 

pip install pandas numpy selenium biopython tqdm
```
3. Install Open Babel using conda or pip. You can also compile Open Babel [from source](https://openbabel.org/docs/dev/Installation/install.html). The goal is to make `obabel` callable from the terminal at least from your current environment.
```
conda install -c conda-forge openbabel
# or
pip install openbabel
```
4. Install `seqfold` package with `pip` and check that `seqfold` is callable from the terminal:
```
pip install seqfold
seqfold
``` 
5. Install `rnafold` within `viennarna` package. It can be done easily with conda. You can also build it [from source](https://github.com/ViennaRNA/ViennaRNA), just make sure that after installation `rnafold` is callable from the terminal.
```
conda install -c bioconda viennarna 
rnafold
```
6. Install `NUPACK`. Source code is available [online](http://www.nupack.org/downloads/source) free of charge after registration. Download the tar bundle and execute the following commands:
```
pip install -U nupack -f path/to/nupack-version/package
```
Additional installation inctructions are available at the [package website](https://docs.nupack.org/start/#maclinux-installation). `NUPACK` is somewhat optional and the pipeline can be completed without it. However the package provides an additional option for secondary structure prediction and screening.

7. Install `fpocket` with conda. Otherwise you can build it [from source](https://github.com/Discngine/fpocket), but make sure that it is callable from the terminal.
```
conda install -c conda-forge fpocket
```
8. `Autodock Vina` is avalibale as a [compiled executable](https://github.com/ccsb-scripps/AutoDock-Vina/releases). Download and place the executable any place on your system that is convinient to you and check that `vina` is active:
```
cd path/to/vina/folder
./vina
```

To start a new workflow, clone this repository locally and rename it as you wish.  
Inspect `export_vars_and_paths.sh` and specify absolute paths to some executables.  
Before running the workflow scripts execute:
```
source export_vars_and_paths.sh
```
This is important for scripts to work properly and find the executables. 

## 1. Run secondary structure prediction

This step is automated with `predict_secondary_structure.py`. There is a help section:
```
cd your_aptafold_job_folder
scripts/predict_secondary_structure.py -h
```

Secondary structure prediction can be done on sequences from fasta file. This is relevant, when you are screening sequencing data from an aptamer selection experiment to a particular target. Provide fasta file, target name and target SMILES notation (can be looked up in PubChem or generated from structure [online](http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html)) and prediction algorithms to use. There is a test fasta file to try:

```
scripts/predict_secondary_structure.py -f test_data/test.fasta -l nivalenol -s 'CC1=CC2C(C(C1=O)O)(C3(C(C(C(C34CO4)O2)O)O)C)CO' --seqfold --rnafold --nupack
```

This will produce `1_secondary_structure/folded_aptamers.csv` table with following columns: 

- "aptamer_name", the name of the aptamer
- "sequence", aptamer sequence
- "ligand_name", a respective ligand name
- "ligand_smiles", SMILES notatation of the ligand molecule
- "seqfold", a colum containing the secondary structure in dot-bracket notation named after the tool which generated it, and the "tool_name_ddg" column with predicted folding energy (in kcal/mol)   

Alternatively you can predict secondary structures of aptamers to different targets by providing the csv table as an input. This table should contain following columns: aptamer_name, sequence, ligand_name, ligand_smiles. Use `test.csv` as an example. Tip: such table can be constructed in excel and exported as csv.
 
```
scripts/predict_secondary_structure.py -c test_data/test.csv --seqfold --rnafold --nupack
```

## 2. Generate aptamer pdb models

Aptamer pdb models are prepared with `generate_aptamer_models.py`. This script has a help section available:
```
scripts/generate_aptamer_models.py -h
```

Provide the csv table produced during the previous step and the secondary structure column to use for 3D-model generation. The script requests pdb structures from [RNAComposer server](https://rnacomposer.cs.put.poznan.pl) through the Chrome browser, which is activated in the background.

```
scripts/generate_aptamer_models.py -i 1_secondary_structure/folded_aptamers.csv -f seqfold
```

Note: It may take a couple of minutes to configure and start the browser. Wait for the download messages in log.

## 3. Docking

When aptamer models are generated you can enter the `2_docking` dir and launch docking runs. The script will convert ligands and aptamers to pdbqt format then search for binding pockets in each aptamer with `fpocket`. For each aptamer the respective ligand will be docked in each found pocket. Inspect `docking.flags` if you want to change some docking parameters for Vina.

To start the docking provide the folder with pdb aptamers and the `folded_aptamers.csv` table.
```
cd 2_docking 
../scripts/run_docking.py -h 

../scripts/run_docking.py -a aptamers_pdb -f scripts/docking.flags -l ../1_secondary_structure/folded_aptamers.csv
```

When the docking finishes, the resulting scores are saved in `docking_data.csv`. The resulting docking positions are available in `docking_out/pdb` dir.
