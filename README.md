# surfaces

clang Vcontacts-v1-2.c -o vcon

### Outdated
For numeric interface interactions between chain ABC... and DE... of pdb file

python surface_cont.py -f file.pdb -c1 ABC -c2 DE -o output.csv

For list of interactions between chains ABC... and DE... of pdb file

python surface_list.py -f file.pdb -c1 ABC -c2 DE -o output.csv

For list of interactions with partial CF values between chains ABC... and DE... of pdb file

python surface_values.py -f file.pdb -c1 ABC -c2 DE -o output.csv

### Up to date

For numeric interface interactions between chain ABC... and DE... of pdb file

python new_surface_cont.py -f file.pdb -c1 ABC -c2 DE -o output.csv -def atomtypes_definition.def -dat atomtypes_interactions.dat

For a visual output as a pymol session

python image_surfaces.py -f pdb_file.pdb -c csv_file.csv -o pymol_session_output.pse
