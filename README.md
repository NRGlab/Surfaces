# surfaces

clang Vcontacts-v1-2.c -o vcon

For numeric interface interactions between chain A and B of pdb file
python surface_cont.py -f file.pdb -c1 A -c2 B -o output.csv

For list of interactions between chains ABC... and DE... of pdb file
python surface_list.py -f file.pdb -c1 ABC -c2 DE -o output.csv
