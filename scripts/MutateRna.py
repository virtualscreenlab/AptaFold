# a scrip for RNA.pdb mutation to DNA.pdb
from Bio import PDB
from Bio.PDB.PDBIO import PDBIO
import numpy as np
from math import sqrt

# build reference from DT.pdb
def get_reference_methyl():
    parser = PDB.PDBParser()
    residue = parser.get_structure('DT', './scripts/DT_reference.pdb')
    ref = {}
    reference_atoms = ["C5", "C7", "H71", "H72", "H73", "H2\'2"]
    for atom in residue.get_atoms():
        if atom.name in reference_atoms:
            ref[atom.name] = atom
    return ref

# retrieve atom in a particular residue
def get_atom(residue, name):
    for atom in residue.get_atoms():
        if atom.name == name:
            return atom
    return None


# a function to calculate H2'2 coordinates from H2', C2' and O2' coordinates
def calculate_H2_coordinates(H2, C2, O2):
    l1 = np.linalg.norm(np.subtract(C2, H2))
    l2 = np.linalg.norm(np.subtract(C2, O2))
    l = l1 / (l2 - l1)
    coordinates = [(C2[0] + l * O2[0]) / (1 + l),
                   (C2[1] + l * O2[1]) / (1 + l),
                   (C2[2] + l * O2[2]) / (1 + l)]
    return coordinates


# a function to calculate coordinates for methyl hydrogens
def calculate_H7_coordinates(C5, C7):
    l = 1.2464789 # precalculated parameters
    r = 1.0364275

    projection = ((C7 - C5) * l) + C5
    p = C7 - C5
    p = p / np.linalg.norm(p)
    p1 = np.array([p[2], 0, -p[0]])
    p1 = p1 / np.linalg.norm(p1)
    p2 = np.cross(p, p1)

    h71 = projection + r * (1 * p1 + 0 * p2)
    h72 = projection + r * (-0.5 * p1 + sqrt(3)/2 * p2)
    h73 = projection + r * (-0.5 * p1 - sqrt(3)/2 * p2)

    return [h71, h72, h73]

# a function to calculate coordinates for C5M atom
def calculate_C7_coordinates(C5, N1, N3):

    p13 = (N1 + N3) / 2
    c7 = C5 - p13
    c7 = c7 / np.linalg.norm(c7)
    c7 = C5 + (c7 * 1.5)
    return c7

# conversion of a residue in chain
def convert_residue(residue, reference, length):
    if residue.resname == 'U':
        residue.resname = 'DT'  # name of the residue is changed
        # change Н5 to С5М
        # use coordinates of other atoms in the ring to construct coordinates for methyl group
        C5 = get_atom(residue, 'C5').coord  # these atoms will be used to calculate coordinates for C5M
        N1 = get_atom(residue, 'N1').coord
        N3 = get_atom(residue, 'N3').coord

        C7 = reference['C7'].copy()  # copy atomic data from reference
        C7.coord = calculate_C7_coordinates(C5, N1, N3)  # assign calculated coordinates to C5M
        residue.add(C7)  # add atom to residue

        residue.detach_child('H5')  # delete H5 atom
        # now let's add methyl hydrogens
        h71, h72, h73 = calculate_H7_coordinates(C5, C7.coord)  # calculate new coordinates

        H71 = reference['H71'].copy()  # copy atomic data from reference
        H72 = reference['H72'].copy()  # copy atomic data from reference
        H73 = reference['H73'].copy()  # copy atomic data from reference

        H71.coord = h71  # assign new coordinates
        H72.coord = h72
        H73.coord = h73

        residue.add(H71)  # add hydrogen atoms
        residue.add(H72)
        residue.add(H73)

    elif residue.resname == 'A':
        residue.resname = 'DA'
    elif residue.resname == 'G':
        residue.resname = 'DG'
    elif residue.resname == 'C':
        residue.resname = 'DC'
    else:
        print(residue.resname)
        # пофиксил тут, добавив ValueError иначе компилятор ругается
        raise ValueError('Unrecognized residue')

    #if residue.id[1] == 1:  # detach extra H atoms at the beginning
    #    residue.detach_child('HTER')
    #if residue.id[1] == length:  # detach extra H atoms at the end
    #    residue.detach_child('HCAP')

    # copy and insert H2'2 atom from reference
    H22 = reference["H2\'2"].copy()  # copy atomic data from reference

    # get coordinates of several atoms to calulate coordinates for new H atom substituting O2'
    C2 = get_atom(residue, 'C2\'').coord
    O2 = get_atom(residue, 'O2\'').coord
    H2 = get_atom(residue, 'H2\'').coord
    H22.coord = calculate_H2_coordinates(H2, C2, O2)  # assign O2' coordinates based on H2' and C2'
    # assign correct residue id to the atom????
    residue.add(H22)  # add atom to the residue
    # detach 2'-OH groups
    residue.detach_child('O2\'')
    residue.detach_child('HO2\'')

    return residue

def convert_pdb(file):
    parser = PDB.PDBParser()
    rna = parser.get_structure('rna', file)
    # generate reference structure to use in further modifications
    reference_DT  = get_reference_methyl()

    length = 0
    for _ in rna.get_residues():  # get the length of the sequence
        length = length + 1

    # for each residue perform mutation to DNA
    for residue in rna.get_residues():
        if residue.resname in ['A', 'U', 'C', 'G']:
            convert_residue(residue, reference_DT, length)
    # save the acquired structure
    io = PDBIO()
    io.set_structure(rna)
    io.save(file)
