#Imports
import argparse
import subprocess
import sys
import os
import re
import pandas as pd
import itertools
import time

#Useful dicts
aa = {'C':'CYS', 'D':'ASP', 'S':'SER', 'Q':'GLN', 'K':'LYS', 'I':'ILE', 'P':'PRO', 'T':'THR', 'F':'PHE', 'N':'ASN', 'G':'GLY', 'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP', 'A':'ALA', 'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET'}


#Input of pdb file and the 2 chains between which we want to evaluate the interactions

def read_residues(pdb_file, chain1, chain2):
    list_chain1 = []
    list_chain2 = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:4] == 'ATOM':
            if line[21] in chain1:
                res_num = re.findall('\d',line[22:27])
                res_num = int(''.join(res_num))
                res_name = line[17:20]
                string = res_name + str(res_num) + line[21]
                if string not in list_chain1:
                    list_chain1.append(string)
            if line[21] in chain2:
                res_num = re.findall('\d',line[22:27])
                res_num = int(''.join(res_num))
                res_name = line[17:20]
                string = res_name + str(res_num) + line[21]
                if string not in list_chain2:
                    list_chain2.append(string)
    chains = chain1 + chain2
    atoms_init = [0]*len(chains)
    atoms_end = [0]*len(chains)
    for line in Lines:
        for i in range(len(chains)):
            if line[:4] == 'ATOM' and line[21] == chains[i]:
                if atoms_init[i] == 0:
                    atoms_init[i] = (int(line[6:12]))
                if int(line[6:12]) > atoms_end[i]:
                    atoms_end[i] = (int(line[6:12]))
    return (list_chain1, list_chain2, chains, atoms_init, atoms_end)

def clean_pdb(pdb_file, chain1, chain2, new_file):
    f = open(pdb_file, 'r')
    Lines1 = f.readlines()
    Lines2 = []
    for line in Lines1:
        if line[:4] == 'ATOM':
            if line[21] in chain1 or line[21] in chain2:
                #print (line)
                Lines2.append(line)
    f.close()
    e = open(new_file, 'w')
    for line in Lines2:
        e.write(line)
    e.close()
    return


#Function to generate the file with the output of perl code vcont.pl

def vcont(pdb_name):
    string = "perl -w vcont.pl -f " + pdb_name + " > vcont_output.txt"
    #print (string)
    os.system(string)
    #print ('vcont done')

def vcon(pdb_name):
    string = "./vcon " + pdb_name + " > vcon_file.txt"
    #print (string)
    os.system(string)
    #print ('vcon done')
    
    
#Functions to fix the names of the chains

def test_chain(atom, chain1, chain2, inits, ends):
    chains = chain1 + chain2
    for i in range(len(chains)):
        if atom >=inits[i] and atom <=ends[i]:
            chain = chains[i]
    return (chain)

def fix_chain(file, chain1, chain2, inits, ends):
    f = open(file, 'r')
    Lines1 = f.readlines()
    Lines2 = []
    for line in Lines1:
        if line[:1] != '#' and line[6:7] == '|':
            chain = test_chain(int(line[:6]), chain1, chain2, inits, ends)
            #print (chain)
            line = line[:31] + chain + line[32:]
            #print (line)
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
    attype = line[10:13]
    resnum = int(line[15:20])
    res = line[25:28]
    chain = line[31:32]
    return (atnum,attype,resnum,res,chain)

def read_surface(line):
    ind = line[35:].index('|')
    surf = (float(line[35:35+ind]))
    return (surf)

def read_interactions(file, matrix, chain1, chain2):
    f = open(file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:1] != '#' and line[6:7] == '|':
            if line[35:38] == 'Sol':
                main_line = line
                main_atnum,main_attype,main_resnum,main_res,main_chain = read_atom(main_line)
                main_residue = main_res+str(main_resnum)+main_chain
            else:
                atnum,attype,resnum,res,chain = read_atom(line)
                other_residue = res+str(resnum)+chain
                surf = read_surface(line)
                if (chain in chain1 and main_chain in chain2) or (chain in chain2 and main_chain in chain1):
                    if main_chain in chain2:
                        #print (other_residue, main_residue)
                        matrix.loc[other_residue, main_residue] += (surf * score(main_attype, main_res, attype, res))/2
                    if main_chain in chain1:
                        #print (main_residue, other_residue)
                        matrix.loc[main_residue, other_residue] += (surf * score(main_attype, main_res, attype, res))/2
  
    return(matrix)

#get the atom type number from file AMINO.def
def atomtype_num(file, res, attyp):
    if attyp == 'OXT':
        attyp = 'O  '
    f = open(file, 'r')
    Lines = f.readlines()
    for i in range(len(Lines)):
        if Lines[i][:10] == 'RESIDU '+res:
            #print(Lines[i][:10], attyp)
            for j in range(1,15):
                #print(Lines[i+j][13:16])
                if Lines[i+j][:6] == 'ATMTYP' and Lines[i+j][13:16] == attyp:
                    attype_num = int(Lines[i+j][10:13])
    return(attype_num)

#get the interaction between atom type 1 and atom type 2 from flexaid file MC_st0r5.2_6.dat
def interactions(file, at1, at2):
    if len(str(at1)) == 1:
        at1 = ' ' + str(at1)
    if len(str(at2)) == 1:
        at2 = ' ' + str(at2)
    f = open(file, 'r')
    Lines = f.readlines()
    #print (str(at1)+'-'+str(at2))
    for line in Lines:
        if str(at1)+'-'+str(at2) == line[5:10] or str(at2)+'-'+str(at1) == line[5:10]:
            interact = float(line[13:])
    f.close()
    return(interact)

#get the final score based on atom type num and interactions
def score(attype1, res1, attype2, res2):
    at1 = atomtype_num('AMINO.def', res1, attype1)
    at2 = atomtype_num('AMINO.def', res2, attype2)
    value = interactions('MC_st0r5.2_6.dat', at1, at2)
    return (value)


#Function to convert sequence in residue names

def seq_to_res(seq, init_posit, dict, chain):
    residues = []
    for res in seq:
        residues.append(dict[res] + str(init_posit) + chain)
        init_posit = init_posit + 1
    return (residues)

#Order interactions and get the list of residues interacting (colunm) for each residue of interest (row)

def list_interactions(index, matrix, res2):
    list_interact = ''
    row = matrix.loc[index , : ]
    row = row.tolist()
    row_sorted = sorted(row, reverse=False)
    for i in row_sorted:
        ind = row.index(i)
        if i < 0:
            list_interact = list_interact + res2[ind] + '(+), '
        if i > 0:
            list_interact = list_interact + res2[ind] + '(-), '
    return (list_interact)



def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-c1","--chain1", action="store")
    parser.add_argument("-c2","--chain2", action="store")
    parser.add_argument("-o","--output_name", action="store")
    args=parser.parse_args()

    #update chain1 and chain2 to be lists of chains instead of only one chain each
    res1, res2, chains, inits, ends = read_residues(args.pdb_file, args.chain1, args.chain2)
    clean_pdb(args.pdb_file, args.chain1, args.chain2, 'clean.pdb')
        
    vcon('clean.pdb')
    fix_chain('vcon_file.txt', args.chain1, args.chain2, inits, ends)
  
    matrix = [ [ 0 for i in range(len(res2)) ] for j in range(len(res1)) ]
    matrix = pd.DataFrame(matrix)
    matrix.columns = res2 
    matrix.index = res1 
        
    matrix = read_interactions('vcon_file.txt', matrix, args.chain1, args.chain2)
    
    #print (matrix)
    
    surface_list = []
    for residue in res1:
        list_interact = list_interactions(residue, matrix, res2)
        surface_list.append(list_interact[:-2])
        if list_interact != '':
            print (residue + ': ' + list_interact[:-2])
    
    matrix['Interactions'] = surface_list
    
    matrix[['Interactions']].to_csv(args.output_name)
    
    
    return
    
main()