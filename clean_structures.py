#Imports
import argparse

def get_atoms(line): #take strings after the ' ' and before the ':'
    list_atoms = []
    L = line[5:-2].split(",")
    for item in L:
        item = item.strip()
        id = item.index(":")
        list_atoms.append(item[:id])
    return (list_atoms)

def read_atom(line):
    atnum = int(line[6:11])
    attype = line[11:17].strip()
    resnum = line[23:30].strip()
    res = line[17:20].strip()
    chain = line[21:22].strip()
    return (atnum,attype,resnum,res,chain)

def check_line(line,res,atoms): #see if the atom of that line is described in the def file
    atnum,attype,resnum,residue,chain = read_atom(line)
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
    
def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-def","--atomtypes_definition", action="store")
    args=parser.parse_args()
    
    res = []
    atoms = []
    
    def_file = open(args.atomtypes_definition, "r")
    Lines = def_file.readlines()
    for line in Lines:
        res.append(line[:3])
        atoms.append(get_atoms(line))
    def_file.close()
    
    #print (res)
    #print (atoms)
    
    pdb = open(args.pdb_file, "r")
    clean_pdb = open("clean_" + args.pdb_file, "w")
    Lines = pdb.readlines()
    for line in Lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            if check_line(line,res,atoms):
                clean_pdb.write(line)
            
    pdb.close()
    clean_pdb.close()
    
    return


main()
