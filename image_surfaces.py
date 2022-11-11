import pymol
import pandas as pd
from colour import Color
import matplotlib.colors

def get_sum_per_residue(surfaces_file):
    residues = []
    values = []
    surf = pd.read_csv(surfaces_file, index_col=0)
    sum_column = surf.sum()
    sum_row = surf.sum(axis='columns')
    for x in range(len(sum_column)-2):
        values.append(sum_column[x])
        residues.append(surf.columns[x])
    for y in range(len(sum_row)-2):
        values.append(sum_row[y])
        residues.append(surf.index[y])
    return (residues, values)

def get_pairs_contacts(surfaces_file):
    pairs = []
    values = []
    surf = pd.read_csv(surfaces_file, index_col=0)
    for i in range (len(surf.index)):
        for j in range (len(surf.columns)):
            if surf.loc[surf.index[i],surf.columns[j]] != 0:
                pairs.append([surf.index[i],surf.columns[j]])
                values.append(surf.loc[surf.index[i],surf.columns[j]])
    return (pairs, values)

def read_residue(res):
    type_res = res[:3]
    chain_res = res[-1]
    num_res = res[3:-1]
    return (type_res, chain_res, num_res)

def color_residue(res, color):
    pymol.cmd.delete('sele')
    type_res, chain_res, num_res = read_residue(res)
    selection_string = 'chain' + chain_res + ' and chain ' + chain_res + ' and resi ' + num_res
    pymol.cmd.set_color(res, color)
    pymol.cmd.select(selection_string)
    #pymol.cmd.show('spheres', 'sele')
    pymol.cmd.set("cartoon_transparency", 0.00, 'sele')
    pymol.cmd.color(res, 'sele')
    return

def color_scale(values):
    
    Total_colors = []
    
    for i in range(5):
        c = Color("red", saturation=1/(i+1))
        Total_colors.append(c.rgb)
    white = Color("white")
    Total_colors.append(white.rgb)
    for i in range(5):
        c = Color("blue", saturation=1/(5-i))
        Total_colors.append(c.rgb)
    print (Total_colors)
    
    max_value = max(values)
    min_value = min(values)
    if abs(min_value) > abs(max_value):
        range_value = 2 * abs(min_value)
    else:
        range_value = 2 * abs(max_value)
    step_value = range_value/10
    
    color_codes = []
    for value in values:
        s = range_value/2 - (-1*value)
        n = int(s // step_value)
        color_codes.append(list(Total_colors[n]))
        
    return (color_codes)

def split_states(residues, pdb_file):
    chains = []
    for res in residues:
        type_res, chain_res, num_res = read_residue(res)
        if chain_res not in chains:
            chains.append(chain_res)
    for C in chains:
        pymol.cmd.select('chain ' + C)
        pymol.cmd.extract('chain' + C, 'sele')
    pymol.cmd.delete(pdb_file[:-4])
    return (chains)

def show_separate_surfaces(chains):
    for C in chains:
        pymol.cmd.show('surface', 'chain' + C)
    return
    
def generate_session(pdb_file, surfaces_file, session_file_name):
    residues, values = get_sum_per_residue(surfaces_file)
    color_codes = color_scale(values)
    pymol.cmd.load(pdb_file)
    pymol.cmd.color('grey60', pdb_file[:-4])
    chains = split_states(residues, pdb_file)
    for C in chains:
        pymol.cmd.set("cartoon_transparency", 0.55, 'chain' + C)
    for i in range(len(residues)):
        if values[i] != 0:
            color_residue(residues[i], color_codes[i])
    show_separate_surfaces(chains)
    pymol.cmd.save(session_file_name, format='pse')
    return

### test function get_sum_per_residue
#residues, values = get_sum_per_residue('LIGIN_output_example.csv')
#print (color_scale(values))

### test function get_pairs_contacts
#pairs, values = get_pairs_contacts('LIGIN_output_example.csv')
#print (pairs)
#print (values)

### test function read_residue
#print (read_residue('TYR513B'))

def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-c","--input_csv_file", action="store")
    parser.add_argument("-o","--pymol_session_output_name", action="store")
    args=parser.parse_args()

    generate_session(args.pdb_file, args.input_csv_file, args.pymol_session_output_name)
    
    return

main()
