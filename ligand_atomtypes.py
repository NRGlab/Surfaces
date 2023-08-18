# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)

import pymol
import pandas as pd
import argparse

# ANY ATOM NOT IN THIS LIST WILL BE CONSIDERED DUMMY:39
dict_types={'C.1':'1', 'C.2':'2', 'C.3':'3', 'C.ar':'4', 'C.cat':'5', 'N.1':'6', 
 'N.2':'7', 'N.3':'8', 'N.4':'9', 'N.ar':'10', 'N.am':'11', 'N.pl3':'12',
 'O.2':'13', 'O.3':'14', 'O.co2':'15', 'O.ar':'16', 'S.2':'17', 'S.3':'18',
 'S.o':'19', 'S.O2':'20', 'S.ar':'21', 'P.3':'22', 'F':'23', 'Cl':'24', 'Br':'25', 'I':'26', 'Se':'27', 'Mg':'28', 'Sr':'29', 
 'Cu':'30', 'Mn':'31', 'Hg':'32', 'Cd':'33', 'Ni':'34', 'Zn':'35', 'Ca':'36', 'Fe':'37', 'Cooh':'38', 'Du':'39'}

def convertto_mol2 (pdb_file):
    pymol.cmd.load(pdb_file)
    pymol.cmd.save(pdb_file[:-4] + '.mol2')
    return

def convertto_pdb (mol_file):
    pymol.cmd.load(mol_file)
    pymol.cmd.save(mol_file[:-5] + '.pdb')
    return

def read_mol2 (mol_file):
    f = open(mol_file, 'r')
    Lines = f.readlines()
    selected_Lines = []
    # select lines after "@<TRIPOS>ATOM" and before "@<TRIPOS>BOND"
    select = False
    for i in range(len(Lines)):
        if Lines[i] == '@<TRIPOS>BOND\n':
            select = False
        if select:
            selected_Lines.append(Lines[i])
        if Lines[i] == '@<TRIPOS>ATOM\n':
            select = True
    # breaks between columns: '\t'
    df = pd.DataFrame(columns=['number','name','coord1','coord2','coord3','type','n','res','k'])
    for line in selected_Lines:
        add = line[:-2].split("\t")
        df.loc[len(df)] = add
    list_atomnames = df['name'].tolist()
    for k in range(len(list_atomnames)):
        list_atomnames[k] = list_atomnames[k].strip()
    list_atomtypes = df['type'].tolist()
    res = df.loc[0,'res']
    return (list_atomnames, list_atomtypes, res)

def atomtypes_to_numbers (list_atomtypes):
    list_atomnumbers = []
    for i in range(len(list_atomtypes)):
        if list_atomtypes[i] in dict_types:
            list_atomnumbers.append(dict_types[list_atomtypes[i]])
        else:
            list_atomnumbers.append(dict_types['Du'])
    return (list_atomnumbers)

def check_atomns (list_atomnames, list_atomtypes):
    for i in range(len(list_atomnames)):
        for j in range(len(list_atomnames)):
            if list_atomnames[i] == list_atomnames[j] and list_atomtypes[i] != list_atomtypes[j]:
                return (False)
    return (True)
   
def check_atom_names (pdb_file):
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            atom_name = line[12:17]
            atom_name = atom_name.strip()
            if len(atom_name) > 3:
                return (False)
    return (True)

def custom_def_file (initial_def_file, list_atomnames, list_atomnumbers, res):
    f = open(initial_def_file, 'r')
    Lines = f.readlines()
    g = open('custom_' + initial_def_file, 'w')
    new_line = '\n' + res[:3] + ' | '
    for i in range(len(list_atomnames)):
        new_line = new_line + list_atomnames[i] + ':' + list_atomnumbers[i] + ', '
    for line in Lines:
        g.write(line)
    g.write(new_line)
    return


def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-fl","--ligand_pdb_file", action="store")
    parser.add_argument("-def","--atomtypes_definition", action="store")
    args=parser.parse_args()
    
    # EVERY ATOM NAME SHOULD HAVE UP TO 3 CHARACTERS IN ORDER TO AVOID MOL2 CONVERSION ISSUES
    if not check_atom_names (args.ligand_pdb_file):
        print ("WARNING: ATOM NAMES LARGER THAN 3 CHARACTERS - POSSIBLE PROBLEM WITH ATOM TYPE READING")
    
    convertto_mol2 (args.ligand_pdb_file)
    list_atomnames, list_atomtypes, res = read_mol2 (args.ligand_pdb_file[:-4] + '.mol2')
    
    # EVERY ATOM FROM THE LIGAND NEEDS TO HAVE A DIFFERENT NAME; EG. CA,CB... 
    if check_atomns (list_atomnames, list_atomtypes):
        list_atomnumbers = atomtypes_to_numbers (list_atomtypes)
        custom_def_file (args.atomtypes_definition, list_atomnames, list_atomnumbers, res)
    else:
        print ("WARNING: ATOMS WITH DIFFERENT ATOM TYPES AND SAME ATOM NAME")
        
    return
    
    # remove files
    os.remove(args.ligand_pdb_file + ".pdb")
    os.remove(args.ligand_pdb_file + ".mol2")
    
main()
