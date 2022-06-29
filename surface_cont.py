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
    atoms_init = [-1]*len(chains)
    atoms_end = [-1]*len(chains)
    for line in Lines:
        for i in range(len(chains)):
            if line[:4] == 'ATOM' and line[21] == chains[i]:
                if atoms_init[i] == -1:
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
        if line[:1] != '#' and line[31:34] != 'Sol' and line != '\n':
            chain = test_chain(int(line[20:30]), chain1, chain2, inits, ends)
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

def read_interactions(file, matrix, chain1, chain2):
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

#Functions related to interaction with solvent

def create_unbound(pdb_file, chain, new_pdb_file):
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    new_Lines = []
    for line in Lines:
        if line[:4] == 'ATOM' and line[21] in chain:
                new_Lines.append(line)
    f.close()
    e = open(new_pdb_file, 'w') 
    for line in new_Lines: 
        e.write(line)          
    e.close()
    return

def solvent_score(attype, res):
    at1 = atomtype_num('AMINO.def', res, attype)
    at2 = 40
    value = interactions('MC_st0r5.2_6.dat', at1, at2)
    return (value)

def read_solvent_surface(line):
    ind = line.index('\n')
    surf = (float(line[50:50+ind]))
    return (surf)

def interactions_solvent(residues_list, file):
    lista = [0]*len(residues_list)
    f = open(file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if len(line) >=6 and line[31:34] == 'Sol':
            surface = read_solvent_surface(line)
            atnum,attype,resnum,res,chain = read_atom(line)
            residue = res+str(resnum)+chain
            if residue in residues_list:
                score = solvent_score(attype, res)
                ind = residues_list.index(residue)
                lista[ind] = lista[ind] + (score*surface)
    return(lista)




def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-c1","--chain1", action="store")
    parser.add_argument("-c2","--chain2", action="store")
    parser.add_argument("-o","--output_name", action="store")
    args=parser.parse_args()

    res1, res2, chains, inits, ends = read_residues(args.pdb_file, args.chain1, args.chain2)
    clean_pdb(args.pdb_file, args.chain1, args.chain2, 'clean.pdb')
    
    #results for the specific mut

    vcon('clean.pdb')
    fix_chain('vcon_file.txt', args.chain1, args.chain2, inits, ends)
  
    matrix = [ [ 0 for i in range(len(res2) + 2) ] for j in range(len(res1) + 2) ]
    matrix = pd.DataFrame(matrix)
    matrix.columns = res2 + ['solvent bound', 'solvent unbound']
    matrix.index = res1 + ['solvent bound', 'solvent unbound']
    
    matrix = read_interactions('vcon_file.txt', matrix, args.chain1, args.chain2)
    
    print (matrix)
    
    
    #generating interactions of RBD with solvent and ACE2 with solvent for bound state
    
    boundRBD = interactions_solvent(res1, 'vcon_file.txt')
    boundACE2 = interactions_solvent(res2, 'vcon_file.txt')
    
    #creating an unbound state of RBD, generating interactions of RBD with solvent for unbound state

    create_unbound('clean.pdb', args.chain1,'unbound1.pdb')
    vcon('unbound1.pdb')
    fix_chain('vcon_file.txt', args.chain1, args.chain2, inits, ends)
    unboundRBD = interactions_solvent(res1, 'vcon_file.txt')
    
    #generating interactions of ACE2 with solvent for unbound state - done only once, result now put as a row that will be added to every mutant
    
    create_unbound(args.pdb_file, args.chain2,'unbound2.pdb')
    vcon('unbound2.pdb')
    fix_chain('vcon_file.txt', args.chain1, args.chain2, inits, ends)
    unboundACE2 = interactions_solvent(res2, 'vcon_file.txt')
    #print (unboundACE2)
    
    #updating the matrix
    
    matrix['solvent bound'] = boundRBD + [0,0]
    matrix['solvent unbound'] = unboundRBD + [0,0]
    matrix.loc['solvent bound',:] = boundACE2 + [0,0]
    matrix.loc['solvent unbound',:] = unboundACE2 + [0,0]
    
    print (matrix)
    
    matrix.to_csv(args.output_name)

    
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
