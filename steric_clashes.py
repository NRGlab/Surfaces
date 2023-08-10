# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)

import argparse
import re


def get_atoms(pdb_file):
    atoms = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            atom_n = (int(line[6:12]))
            coord1 = (float(line[31:38]))
            coord2 = (float(line[39:46]))
            coord3 = (float(line[47:54]))
            #print (coord1, coord2, coord3)
            atoms.append([atom_n, coord1, coord2, coord3])
    return (atoms)

def get_residue_name(pdb_line):
    res_num = re.findall('\d',pdb_line[22:27])
    res_num = int(''.join(res_num))
    res_name = pdb_line[17:20]
    string = res_name + str(res_num) + pdb_line[21]
    return (string, res_num)

def get_residues(pdb_file):
    residue_names = []
    residue_numbers = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            name, res_num = get_residue_name(line)
            residue_names.append(name)
            residue_numbers.append(res_num)
    return (residue_names, residue_numbers)

def check_bond(atom1, atom2):
    pymol.cmd.select('id ' + str(atom1))
    pymol.cmd.set_name('sele', 'atom1')
    pymol.cmd.select('bound_to atom1')
    list_bound = pymol.cmd.identify('sele')
    #print (list_bound)
    if atom2 in list_bound:
        return (True)
    else:
        return (False)

def calculate_distance(atom1, atom2):
    delta1 = abs(atom1[1] - atom2[1])
    delta2 = abs(atom1[2] - atom2[2])
    delta3 = abs(atom1[3] - atom2[3])
    dist = ((delta1**2) + (delta2**2) + (delta3**2))**(0.5)
    #print (dist)
    if dist > 2:
        return (False)
    else:
        return (True)
    
def check_adjacence(res_num1, res_num2):
    if abs(res_num1 - res_num2) > 1:
        return (False)
    else:
        return (True)

def give_warning(pairs):
    
    for pair in pairs:
        print ("WARNING: STERIC CLASH BETWEEN RESIDUES " + pair[0] + " AND " + pair[1])
    
    return

def get_steric_clashes(pdb_file):
    
    # get numbers of atoms
    atoms = get_atoms(pdb_file)
    # get names and numbers of residues
    residue_names, residue_numbers = get_residues(pdb_file)

    
    # iterate over every atom
    pairs = []
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            if j > i:
                if (calculate_distance(atoms[i], atoms[j])) and not (check_adjacence(residue_numbers[i], residue_numbers[j])):
                    pair = [residue_names[i], residue_names[j]]
                    if pair not in pairs:
                        pairs.append(pair)
    
    give_warning (pairs)

    return

def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    args=parser.parse_args()
    
    get_steric_clashes(args.pdb_file)
    
    return

main()
