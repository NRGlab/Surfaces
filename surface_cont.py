# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)

#Imports
import argparse
import sys
import os
import re
import pandas as pd
from vconpy import run_vcon


#Useful dicts
aa = {'C':'CYS', 'D':'ASP', 'S':'SER', 'Q':'GLN', 'K':'LYS', 'I':'ILE', 'P':'PRO', 'T':'THR', 'F':'PHE', 'N':'ASN', 'G':'GLY', 'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP', 'A':'ALA', 'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET'}

#Input of pdb file and the 2 chains between which we want to evaluate the interactions

def read_residues(pdb_file, chain1, chain2):
    list_chain1 = []
    list_chain2 = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            if line[21] in chain1:
                res_num = re.findall(r'[+-]?\d+',line[22:27])
                res_name = line[17:20]
                string = res_name + str(res_num[0]) + line[21]
                if string not in list_chain1:
                    list_chain1.append(string)
            if line[21] in chain2:
                res_num = re.findall(r'[+-]?\d+',line[22:27])
                res_name = line[17:20]
                string = res_name + str(res_num[0]) + line[21]
                if string not in list_chain2:
                    list_chain2.append(string)
    chains = chain1 + chain2
    atoms_numbers = []
    for k in range(len(chains)):
        atoms_numbers.append([])
    for line in Lines:
        for i in range(len(chains)):
            if (line[:4] == 'ATOM' or line[:4] == 'HETA') and line[21] == chains[i]:
                atoms_numbers[i].append(int(line[6:12]))
    return (list_chain1, list_chain2, chains, atoms_numbers)

#Function to generate the file with the output of perl code vcont.pl

def vcon(pdb_name):
    #run_vcon(pdb_name, showbonded=True)
    os.system(f"./vcon {pdb_name} > {pdb_name[:-4]}.vcon")

def get_chain(atom, og_chain, chains, atom_numbers):
    if og_chain in chains:
        index = chains.index(og_chain)
    else:
        return (og_chain)
    if atom in atom_numbers[index]:
        return (og_chain)
    else:
        for i in range(len(chains)):
            if i != index:
                if atom in atom_numbers[i]:
                    return (chains[i])
    return (0)

#Functions to open vcon results, .def and .dat for atom types interactions

def read_atom(line):
    atnum = int(line[:6])
    attype = line[8:11]
    resnum = int(line[12:17])
    res = line[19:22]
    chain = line[23:24]
    return (atnum,attype,resnum,res,chain)

def read_surface(line):
    # Extract all floats on the line (handles 0, 0., 1.23, 1e-3, etc.)
    floats = re.findall(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?', line)
    if 'Sol_acc_surf' in line:
        # Atom header line: only one numeric value at the end (the solvent-accessible area)
        return float(floats[-1])
    else:
        # Contact line: last two numbers are AREA and DIST → take penultimate as AREA
        if len(floats) < 2:
            raise ValueError(f"Could not find AREA/DIST on line: {line!r}")
        return float(floats[-2])

def read_interactions(file, matrix, chain1, chain2, def_file, dat_file, atom_numbers, scale_factor):
    chains = chain1 + chain2
    f = open(file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:1] != '#' and line != '\n':
            if line[31:34] == 'Sol':
                main_line = line
                main_atnum,main_attype,main_resnum,main_res,main_chain = read_atom(main_line)
                fixed_main_chain = get_chain(main_atnum,main_chain,chains,atom_numbers)
                if fixed_main_chain != 0:
                    main_residue = main_res+str(main_resnum)+fixed_main_chain
            else:
                atnum,attype,resnum,res,chain = read_atom(line[22:])
                fixed_other_chain = get_chain(atnum,chain,chains,atom_numbers)
                if fixed_other_chain != 0:
                    other_residue = res+str(resnum)+fixed_other_chain
                    surf = read_surface(line)
                    if (fixed_other_chain in chain1 and fixed_main_chain in chain2) or (fixed_other_chain in chain2 and fixed_main_chain in chain1):
                        if fixed_main_chain in chain2:
                            matrix.loc[other_residue, main_residue] += (surf * score(main_attype, main_res, attype, res, def_file, dat_file) * scale_factor)/2
                        if fixed_main_chain in chain1:
                            matrix.loc[main_residue, other_residue] += (surf * score(main_attype, main_res, attype, res, def_file, dat_file) * scale_factor)/2
    return(matrix)

#get the atom type number from def file of choice
def atomtype_num(def_file, res, attyp):
    attyp = attyp.replace(" ", "")
    f = open(def_file, 'r')
    Lines = f.readlines()
    for i in range(len(Lines)):
        if Lines[i][:3] == res:
            ind = Lines[i].index(attyp+':')
            ind_end = Lines[i][ind:].index(',')
            attype_num = int(Lines[i][ind+len(attyp)+1:ind+ind_end])
    f.close()
    return(attype_num)

#get the interaction between atom type 1 and atom type 2 from dat file of choice
def interactions(dat_file, at1, at2):
    if len(str(at1)) == 1:
        at1 = ' ' + str(at1)
    if len(str(at2)) == 1:
        at2 = ' ' + str(at2)
    f = open(dat_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if str(at1)+'-'+str(at2) == line[5:10] or str(at2)+'-'+str(at1) == line[5:10]:
            interact = float(line[13:])
    f.close()
    return(interact)

#get the final score based on atom type num and interactions
def score(attype1, res1, attype2, res2, def_file, dat_file):
    at1 = atomtype_num(def_file, res1, attype1)
    at2 = atomtype_num(def_file, res2, attype2)
    value = interactions(dat_file, at1, at2)
    return (value)

#create file of list of interactions
def list_file(matrix,output_name):
    residues1 = []
    residues2 = []
    values = []
    abs_values = []
    for i in range(len(matrix.index)):
        for j in range(len(matrix.columns)):
            num = matrix.loc[matrix.index[i], matrix.columns[j]]
            if num != 0:
                residues1.append(matrix.index[i])
                residues2.append(matrix.columns[j])
                values.append(num)
                abs_values.append(abs(num))
    sorted_residues1 = [x for _,x in sorted(zip(abs_values,residues1),reverse=True)]
    sorted_residues2 = [x for _,x in sorted(zip(abs_values,residues2),reverse=True)]
    sorted_values = [x for _,x in sorted(zip(abs_values,values),reverse=True)]
    f = open(output_name[:-4] + "_list.txt", "w")
    for k in range(len(values)):
        f.write(sorted_residues1[k] + "," + sorted_residues2[k] + "," + str(sorted_values[k]) + "\n")
    f.close()
    return

def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-c1","--chain1", action="store")
    parser.add_argument("-c2","--chain2", action="store")
    parser.add_argument("-o","--output_name", action="store")
    parser.add_argument("-def","--atomtypes_definition", action="store")
    parser.add_argument("-dat","--atomtypes_interactions", action="store")
    args=parser.parse_args()
    
    res1, res2, chains, atom_numbers = read_residues(args.pdb_file, args.chain1, args.chain2)
    #print (res1, res2)
    #print (chains, atom_numbers)
        
    vcon(args.pdb_file)
  
    matrix = [ [ 0.0 for i in range(len(res2)) ] for j in range(len(res1)) ]
    matrix = pd.DataFrame(matrix)
    matrix.columns = res2
    matrix.index = res1
    
    # Determined according to the AB-Bind dataset results
    scale_factor = 0.00024329
    
    matrix = read_interactions(args.pdb_file[:-4]+".vcon", matrix, args.chain1, args.chain2, args.atomtypes_definition, args.atomtypes_interactions, atom_numbers, scale_factor)

    matrix.to_csv(args.output_name)

    list_file(matrix,args.output_name)
    
    os.remove(args.pdb_file[:-4]+".vcon")
    
    return
    
main()
