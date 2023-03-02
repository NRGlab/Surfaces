import argparse

def read_atom(line):
    atnum = int(line[6:11])
    attype = line[11:16].strip()
    resnum = line[23:30].strip()
    res = line[17:20].strip()
    chain = line[21:22].strip()
    return (atnum,attype,resnum,res,chain)

def create_ligand_file (pdb_file_name, lig_name):
    pdb_file = open(pdb_file_name, "r")
    lig_pdb_file = open(lig_name+".pdb", "w")
    Lines = pdb_file.readlines()
    for line in Lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            atnum,attype,resnum,res,chain = read_atom(line)
            if res == lig_name:
                lig_pdb_file.write(line)
    lig_pdb_file.close()
    return

def main():

    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-lig","--lig_name", action="store")
    args=parser.parse_args()
    
    create_ligand_file (args.pdb_file, args.lig_name)

    return
    
main()
