import argparse
import pandas as pd
import sys
import os
import re

def read_atom(line):
    atnum = int(line[6:11])
    attype = line[11:16].strip()
    resnum = line[23:30].strip()
    res = line[17:20].strip()
    chain = line[21:22].strip()
    return (atnum,attype,resnum,res,chain)

def create_ligand_file (pdb_file_name, lig_name, lig_out_file):
    pdb_file = open(pdb_file_name, "r")
    lig_pdb_file = open(lig_out_file, "w")
    Lines = pdb_file.readlines()
    for line in Lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            atnum,attype,resnum,res,chain = read_atom(line)
            if res == lig_name:
                lig_pdb_file.write(line)
    lig_pdb_file.close()
    return

# ANY ATOM NOT IN THIS LIST WILL BE CONSIDERED DUMMY:39
dict_types={'C.1':'1', 'C.2':'2', 'C.3':'3', 'C.ar':'4', 'C.cat':'5', 'N.1':'6', 
 'N.2':'7', 'N.3':'8', 'N.4':'9', 'N.ar':'10', 'N.am':'11', 'N.pl3':'12',
 'O.2':'13', 'O.3':'14', 'O.co2':'15', 'O.ar':'16', 'S.2':'17', 'S.3':'18',
 'S.o':'19', 'S.O2':'20', 'S.ar':'21', 'P.3':'22', 'F':'23', 'Cl':'24', 'Br':'25', 'I':'26', 'Se':'27', 'Mg':'28', 'Sr':'29', 
 'Cu':'30', 'Mn':'31', 'Hg':'32', 'Cd':'33', 'Ni':'34', 'Zn':'35', 'Ca':'36', 'Fe':'37', 'Cooh':'38', 'Du':'39'}

def convertto_mol2 (pdb_file):
    import pymol
    pymol.cmd.load(pdb_file)
    pymol.cmd.save(pdb_file[:-4] + '.mol2')
    return

def convertto_pdb (mol_file):
    import pymol
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

def check_atoms (list_atomnames, list_atomtypes):
    for i in range(len(list_atomnames)):
        for j in range(len(list_atomnames)):
            if list_atomnames[i] == list_atomnames[j] and list_atomtypes[i] != list_atomtypes[j]:
                return (False)
    return (True)
   
def check_atom_names_pdb (pdb_file):
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            atom_name = line[12:17]
            atom_name = atom_name.strip()
            if len(atom_name) > 3:
                return (False)
    return (True)
    
def check_atom_names_mol2 (mol2_file):
    f = open(mol2_file, 'r')
    Lines = f.readlines()
    selected_Lines = []
    select = False
    for i in range(len(Lines)):
        if Lines[i] == '@<TRIPOS>BOND\n':
            select = False
        if select:
            selected_Lines.append(Lines[i])
        if Lines[i] == '@<TRIPOS>ATOM\n':
            select = True
    df = pd.DataFrame(columns=['number','name','coord1','coord2','coord3','type','n','res','k'])
    for line in selected_Lines:
        add = line[:-2].split("\t")
        df.loc[len(df)] = add
    list_atomnames = df['name'].tolist()
    for atom_name in list_atomnames:
        atom_name = atom_name.strip()
        if len(atom_name) > 3:
            return (False)
    return (True)

def custom_def_file(initial_def_file, list_atomnames, list_atomnumbers, res, output_file_name):
    f = open(initial_def_file, 'r')
    Lines = f.readlines()
    g = open(f"{output_file_name}", 'w')
    new_line = '\n' + res[:3] + ' | '
    for i in range(len(list_atomnames)):
        new_line = new_line + list_atomnames[i] + ':' + list_atomnumbers[i] + ', '
    for line in Lines:
        g.write(line)
    g.write(new_line)
    return

def remove_duplicates(list_atomnames, list_atomnumbers):
    new_list_atomnames = []
    new_list_atomnumbers = []
    for k in range(len(list_atomnames)):
        if list_atomnames[k] not in new_list_atomnames:
            new_list_atomnames.append(list_atomnames[k])
            new_list_atomnumbers.append(list_atomnumbers[k])
    return (new_list_atomnames, new_list_atomnumbers)

def get_atoms(line): #take strings after the ' ' and before the ':'
    list_atoms = []
    L = line[5:-2].split(",")
    for item in L:
        item = item.strip()
        id = item.index(":")
        list_atoms.append(item[:id])
    return (list_atoms)

def read_atom2(line):
    atnum = int(line[6:11])
    attype = line[11:16].strip()
    resnum = line[23:30].strip()
    res = line[17:20].strip()
    chain = line[21:22].strip()
    return (atnum,attype,resnum,res,chain)

def check_line(line,res,atoms): #see if the atom of that line is described in the def file
    atnum,attype,resnum,residue,chain = read_atom2(line)
    if residue in res:
        id = res.index(residue)
        if attype in atoms[id]:
            return (True)
        else:
            print ("WARNING: ATOM " + attype + " NOT DEFINED FOR " + residue)
            return (False)
    else:
        print ("WARNING: " + residue + " NOT DEFINED")
        return (False)
    
#Useful dicts
aa = {'C':'CYS', 'D':'ASP', 'S':'SER', 'Q':'GLN', 'K':'LYS', 'I':'ILE', 'P':'PRO', 'T':'THR', 'F':'PHE', 'N':'ASN', 'G':'GLY', 'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP', 'A':'ALA', 'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET'}

#Input of pdb file and the 2 chains between which we want to evaluate the interactions

def read_residues(pdb_file, chains, ligand):
    list_chains = []
    list_atoms = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            if line[21] in chains:
                res_num = re.findall('[+-]?\d+',line[22:27])
                res_name = line[17:20]
                string = res_name + str(res_num[0]) + line[21]
                if string not in list_chains and res_name != ligand:
                    list_chains.append(string)
            if line[17:20] == ligand:
                atom_name = line[12:16].strip()
                atom_number = re.findall('[+-]?\d+',line[6:12])
                string = str(atom_number[0]) + atom_name
                list_atoms.append(string)
    atoms_numbers = []
    for k in range(len(chains)):
        atoms_numbers.append([])
    for line in Lines:
        for i in range(len(chains)):
            if (line[:4] == 'ATOM' or line[:4] == 'HETA') and line[21] == chains[i]:
                atoms_numbers[i].append(int(line[6:12]))
    return (list_chains, list_atoms, atoms_numbers)

#Function to generate the file with the output of perl code vcont.pl

def vcon(pdb_name, file_out):
    string = f'{os.path.join(".", "vcon")} {pdb_name} > {file_out}'
    os.system(string)

#Functions to fix the names of the chains

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

#Functions to open vcon results, AMINO.def and MC_st0r5.2_6.dat for atom types interactions

def read_atom3(line):
    atnum = int(line[:6])
    attype = line[7:11].strip()
    resnum = int(line[12:17])
    res = line[19:22]
    chain = line[23:24]
    return (atnum,attype,resnum,res,chain)

def read_surface(line):
    surf = (float(line[-12:-6]))
    return (surf)

def read_interactions(file, matrix, chains, ligand, def_file, dat_file, atom_numbers, scale_factor):
    f = open(file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:1] != '#' and line != '\n':
            if line[31:34] == 'Sol':
                atnum,attype,resnum,res,chain = read_atom3(line)
                fixed_chain = get_chain(atnum,chain,chains,atom_numbers)
                if res == ligand:
                    atom_name = str(atnum) + attype
            else:
                main_atnum,main_attype,main_resnum,main_res,main_chain = read_atom3(line[22:])
                fixed_main_chain = get_chain(main_atnum,main_chain,chains,atom_numbers)
                if fixed_main_chain != 0:
                    if main_res != ligand:
                        main_residue = main_res+str(main_resnum)+fixed_main_chain
                        surf = read_surface(line)
                        if (fixed_main_chain in chains) and (res == ligand):
                            matrix.loc[main_residue, atom_name] += (surf * score(main_attype, main_res, attype, res, def_file, dat_file) * scale_factor)
  
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
    f = open(output_name, "w")
    for k in range(len(values)):
        f.write(sorted_residues1[k] + "," + sorted_residues2[k] + "," + str(sorted_values[k]) + "\n")
    f.close()
    return

def main():

    parser= argparse.ArgumentParser(description="the arguments.")
    parser.add_argument("-lig_id","--ligand_id", action="store")
    parser.add_argument("-lig_file","--ligand_file", action="store")
    parser.add_argument("-dir","--pdb_dir", action="store")
    parser.add_argument("-def","--atomtypes_definition", action="store")
    parser.add_argument("-prefix","--def_prefix", action="store", default='custom', help = "prefix for the output def file")
    parser.add_argument("-o","--output_name", action="store")
    parser.add_argument("-outdir","--output_dir", action="store")
    parser.add_argument("-dat","--atomtypes_interactions", action="store")
    parser.add_argument("-vcon","--vcon_out", action="store")
    args=parser.parse_args()

    lig_id = args.ligand_id
    lig_file = pd.read_csv(args.ligand_file)
    pdb_id = lig_file.loc[lig_file.sample_name == lig_id, 'pdb_id'].values[0]
    pdb_file_name = f"{args.pdb_dir}/{pdb_id}.pdb"

    chains = lig_file.loc[lig_file.sample_name == lig_id, 'cath_chain'].values[0]
    lig_name = lig_file.loc[lig_file.sample_name == lig_id, 'bound_ligand_name'].values[0]
    lig_file_out = f"{args.output_dir}/{lig_name}.pdb"

    create_ligand_file(pdb_file_name, lig_name, lig_file_out)
    lig_file_out_mol2 = f"{args.output_dir}/{lig_name}.mol2"
    convertto_mol2(lig_file_out)

    list_atomnames, list_atomtypes, res = read_mol2 (lig_file_out_mol2)
    # EVERY ATOM NAME SHOULD HAVE UP TO 3 CHARACTERS IN ORDER TO AVOID MOL2 CONVERSION ISSUES
    if not check_atom_names_pdb (lig_file_out):
        print ("WARNING: ATOM NAMES LARGER THAN 3 CHARACTERS - POSSIBLE PROBLEM WITH ATOM TYPE READING")
    # EVERY ATOM FROM THE LIGAND NEEDS TO HAVE A DIFFERENT NAME; EG. CA,CB... 
    if not check_atoms (list_atomnames, list_atomtypes):
        print ("WARNING: ATOMS WITH DIFFERENT ATOM TYPES AND SAME ATOM NAME")
    list_atomnumbers = atomtypes_to_numbers (list_atomtypes)
    list_atomnames, list_atomnumbers = remove_duplicates (list_atomnames, list_atomnumbers)
    
    custom_def_file_name = args.output_dir + "/" + args.def_prefix + "_" + args.atomtypes_definition

    custom_def_file(args.atomtypes_definition, list_atomnames, list_atomnumbers, res, custom_def_file_name)
 
    
  
    res = []
    atoms = []
    
    def_file = open(custom_def_file_name, "r")
    Lines = def_file.readlines()
    for line in Lines:
        res.append(line[:3])
        atoms.append(get_atoms(line))
    def_file.close()
    
    #print (res)
    #print (atoms)
    
    pdb_file = open(pdb_file_name, "r")

    clean_pdb_file_name = f"{args.output_dir}/clean_{pdb_id}.pdb"
    clean_pdb_file = open(clean_pdb_file_name, "w")
    Lines = pdb_file.readlines()
    for line in Lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            if check_line(line,res,atoms):
                clean_pdb_file.write(line)
            
    pdb_file.close()
    clean_pdb_file.close()

    res, atoms, atom_numbers = read_residues(clean_pdb_file_name, chains, lig_name)
    #print (res, atoms)
    #print (chains, atom_numbers)
    vcon_out_file = f"{args.output_dir}/{args.vcon_out}"
    vcon(clean_pdb_file_name, vcon_out_file)
  
    matrix = [ [ 0 for i in range(len(atoms)) ] for j in range(len(res)) ]
    matrix = pd.DataFrame(matrix)
    matrix.columns = atoms
    matrix.index = res
    
    # Determined according to the AB-Bind dataset results
    scale_factor = 0.00024329
    
    matrix = read_interactions(vcon_out_file, matrix, chains, lig_name, custom_def_file_name, args.atomtypes_interactions, atom_numbers, scale_factor)
   
    output_file_name = f"{args.output_dir}/{args.output_name}.csv"
    matrix.to_csv(output_file_name)
    
    list_output_name = f"{args.output_dir}/List_{args.output_name}.txt"
    list_file(matrix,list_output_name)
    
    # remove files
    os.remove(vcon_out_file)
    
    return ()
    
main()
