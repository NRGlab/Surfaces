#Imports
import argparse
import subprocess
import sys
import os
import re
import pandas as pd
import itertools
import time
import steric_clashes

#Useful dicts
aa = {'C':'CYS', 'D':'ASP', 'S':'SER', 'Q':'GLN', 'K':'LYS', 'I':'ILE', 'P':'PRO', 'T':'THR', 'F':'PHE', 'N':'ASN', 'G':'GLY', 'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP', 'A':'ALA', 'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET'}

#Input of pdb file and list of residues for us to evaluate the interactions

def get_residue_name(pdb_line):
    res_num = re.findall('\d',pdb_line[22:27])
    res_num = int(''.join(res_num))
    res_name = pdb_line[17:20]
    string = res_name + str(res_num) + pdb_line[21]
    return (string)

def read_chains(pdb_file):
    chains = []
    atoms = []
    residues = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:4] == 'ATOM':
            # correspondent chain for each atom
            chains.append(line[21])
            atoms.append(int(line[6:12]))
            # correct names for each residue    
            string = get_residue_name(line)
            if string not in residues:
                    residues.append(string)              
    return (residues, chains, atoms)

def read_residues(pdb_file, list_residues):
    atoms_numbers = []
    for k in range(len(list_residues)):
        atoms_numbers.append([])
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        for i in range(len(list_residues)):
            if (line[:4] == 'ATOM' or line[:4] == 'HETA') and (get_residue_name(line)==list_residues[i]):
                atoms_numbers[i].append(int(line[6:12]))
    return (atoms_numbers)


#Function to generate the file with the output of vcon

def vcon(pdb_name):
    string = "./vcon " + pdb_name + " > vcon_file.txt"
    #print (string)
    os.system(string)
    #print ('vcon done')


#Functions to fix the names of the chains

def test_chain(atom, chains, atoms):
    for i in range(len(chains)):
        if atom == atoms[i]:
            chain = chains[i]
    return (chain)

def fix_chain(file, chains, atoms):
    f = open(file, 'r')
    Lines1 = f.readlines()
    Lines2 = []
    for line in Lines1:
        if line[:1] != '#' and line[31:34] != 'Sol' and line != '\n':
            chain = test_chain(int(line[20:30]), chains, atoms)
            line = line[:45] + chain + line[46:]
            Lines2.append(line)
        else:
            Lines2.append(line)
    f.close()
    e = open(file, 'w')
    for line in Lines2:
        e.write(line)
    e.close()
    return


#Functions to open vcon results, AMINO.def and MC_st0r5.2_6.dat for atom types interactions

def read_atom(line):
    atnum = int(line[:6])
    attype = line[8:11]
    resnum = int(line[12:17])
    res = line[19:22]
    chain = line[23:24]
    #print (atnum)
    #print (attype)
    #print (resnum)
    #print (res)
    #print (chain)
    return (atnum,attype,resnum,res,chain)

def read_surface(line):
    surf = (float(line[-6:-1]))
    #print (surf)
    return (surf)

def read_interactions(file, matrix, sele_res, all_res, def_file, dat_file):
    f = open(file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:1] != '#' and line != '\n':
            if line[31:34] == 'Sol':
                #print (line)
                main_line = line
                main_atnum,main_attype,main_resnum,main_res,main_chain = read_atom(main_line)
                main_residue = main_res+str(main_resnum)+main_chain
            else:
                #print (line[22:])
                atnum,attype,resnum,res,chain = read_atom(line[22:])
                other_residue = res+str(resnum)+chain
                surf = read_surface(line)
                if main_residue in sele_res:
                    #print (main_residue, other_residue)
                    matrix.loc[main_residue, other_residue] += (surf * score(main_attype, main_res, attype, res, def_file, dat_file))
  
    return(matrix)


#get the atom type number from def file of choice
def atomtype_num(def_file, res, attyp):
    # solve 'OXT' issue
    attyp = attyp.replace(" ", "")
    f = open(def_file, 'r')
    Lines = f.readlines()
    for i in range(len(Lines)):
        if Lines[i][:3] == res:
            print (Lines[i])
            ind = Lines[i].index(attyp+':')
            print (Lines[i][ind:])
            ind_end = Lines[i][ind:].index(',')
            print (Lines[i][ind+len(attyp)+1:ind+ind_end])
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
    #print (str(at1)+'-'+str(at2))
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

#Function to convert sequence in residue names

def seq_to_res(seq, init_posit, dict, chain):
    residues = []
    for res in seq:
        residues.append(dict[res] + str(init_posit) + chain)
        init_posit = init_posit + 1
    return (residues)


def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-res","--list_residues", action="store")
    parser.add_argument("-o","--output_name", action="store")
    parser.add_argument("-def","--atomtypes_definition", action="store")
    parser.add_argument("-dat","--atomtypes_interactions", action="store")
    args=parser.parse_args()
    
    steric_clashes.get_steric_clashes(args.pdb_file)

    sele_res = list(args.list_residues.split(","))
    print (sele_res)
    all_res, chains, atoms = read_chains(args.pdb_file)
    atoms_numbers = read_residues(args.pdb_file, sele_res)
        
    vcon(args.pdb_file)
    fix_chain('vcon_file.txt',chains, atoms)
  
    matrix = [ [ 0 for i in range(len(all_res)) ] for j in range(len(sele_res)) ]
    matrix = pd.DataFrame(matrix)
    matrix.columns = all_res
    matrix.index = sele_res
    
    print(matrix)
    
    matrix = read_interactions('vcon_file.txt', matrix, sele_res, all_res, args.atomtypes_definition, args.atomtypes_interactions)
    
    print (matrix)
    
    matrix.to_csv(args.output_name)
    
    os.remove('vcon_file.txt')
    
    #TESTS
    
    #####TEST FIX CHAIN FUNCION#####
    #fix_chain('test_file.txt')
    
    #####TEST READ ATOM FUNCTION#####
    #f = open('test_file.txt', 'r')
    #Lines = f.readlines()
    #for line in Lines:
    #    if line[:1] != '#' and line[6:7] == '|':
    #        print(read_atom(line))
              
    #####TEST INTERACTIONS TAKEN FROM FLEXAID FILE#####
    #interact1 = interactions('MC_st0r5.2_6.dat', 7, 14)
    #interact2 = interactions('MC_st0r5.2_6.dat', 7, 18)
    #print(interact1*interact2)
    
    #####TEST ATOM TYPE NUMBERS TAKEN FROM FLEXAID FILE#####
    #print (atomtype_num('AMINO.def', 'THR', 'OG1'))
    
    #####TEST ITEMS ON THE MATRIX#####
    #print(matrix.loc['PHE338F', 'THR27D'])
    #matrix.loc['PHE338F', 'THR27D'] += 1
    #print (matrix)
    
    return
    
main()
