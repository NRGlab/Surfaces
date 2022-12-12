import pymol
import re


def get_atoms(pdb_file):
    atom_numbers = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            atom_n = (int(line[6:12]))
            atom_numbers.append(atom_n)
    return (atom_numbers)

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

def check_distance(atom1, atom2):
    pymol.cmd.select('id ' + str(atom1))
    pymol.cmd.set_name('sele', 'atom1')
    pymol.cmd.select('id ' + str(atom2))
    pymol.cmd.set_name('sele', 'atom2')
    dist = pymol.cmd.get_distance('atom1','atom2')
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
    atom_numbers = get_atoms(pdb_file)
    # get names and numbers of residues
    residue_names, residue_numbers = get_residues(pdb_file)
    
    # load file into pymol
    pymol.cmd.load(pdb_file)
    
    # iterate over every atom
    pairs = []
    for i in range(len(atom_numbers)):
        for j in range(len(atom_numbers)):
            if j > i:
                if (check_distance(atom_numbers[i], atom_numbers[j])) and not (check_adjacence(residue_numbers[i], residue_numbers[j])):
                    pair = [residue_names[i], residue_names[j]]
                    if pair not in pairs:
                        pairs.append(pair)
    
    give_warning (pairs)

    return