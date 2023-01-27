#Imports
import argparse
import subprocess
import sys
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
    type_res, chain_res, num_res = read_residue(res)
    selection_string = 'chain' + chain_res + ' and resi ' + num_res
    pymol.cmd.set_color(res, color)
    pymol.cmd.select(selection_string)
    #pymol.cmd.show('spheres', 'sele')
    pymol.cmd.set("cartoon_transparency", 0.00, 'sele')
    pymol.cmd.color(res, 'sele')
    pymol.cmd.delete('sele')
    return

def color_scale(values, color_scale_range):
    
    Total_colors = []
    
    for i in range(5):
        c = Color("blue", saturation=1/(i+1))
        Total_colors.append(c.rgb)
    white = Color("white")
    Total_colors.append(white.rgb)
    for i in range(5):
        c = Color("red", saturation=1/(5-i))
        Total_colors.append(c.rgb)
    #print (Total_colors)
    
    if color_scale_range is None:
        max_value = max(values)
        min_value = min(values)
        if abs(min_value) > abs(max_value):
            range_value = 2 * abs(min_value)
        else:
            range_value = 2 * abs(max_value)
        step_value = range_value/10
    else:
        color_scale_range = list(color_scale_range[1:-1].split(","))
        min_value = float(color_scale_range[0])
        max_value = float(color_scale_range[1])
        range_value = max_value - min_value
        step_value = range_value/10
    
    color_codes = []
    for value in values:
        s = range_value/2 - (-1*value)
        n = int(s // step_value)
        if n < 0:
            n = 0
        elif n > len(Total_colors):
            n = len(Total_colors)
        color_codes.append(list(Total_colors[n]))
        
    return (color_codes)

def color_distance(pair, value, color, selected_pairs):
    #create distance object
    distance_string = 'dashed_' + pair[0] + '-' + pair[1]
    type_res1, chain_res1, num_res1 = read_residue(pair[0])
    type_res2, chain_res2, num_res2 = read_residue(pair[1])
    selection_string = 'chain' + chain_res1 + ' and resi ' + num_res1 + ' and n. CA'
    pymol.cmd.select(selection_string)
    pymol.cmd.set_name('sele', 'res1')
    selection_string = 'chain' + chain_res2 + ' and resi ' + num_res2 + ' and n. CA'
    pymol.cmd.select(selection_string)
    pymol.cmd.set_name('sele', 'res2')
    pymol.cmd.set_color(distance_string, color)
    pymol.cmd.distance(distance_string, 'res1', 'res2')
    pymol.cmd.color(distance_string, distance_string)
    pymol.cmd.hide('labels', distance_string)
    pymol.cmd.delete('res1')
    pymol.cmd.delete('res2')
    if pair not in selected_pairs:
        pymol.cmd.disable(distance_string)
    return
    
def label_pairs(pair,selected_pairs):
    #create selection
    pair_string = pair[0] + '-' + pair[1]
    type_res1, chain_res1, num_res1 = read_residue(pair[0])
    type_res2, chain_res2, num_res2 = read_residue(pair[1])
    selection_string1 = 'chain' + chain_res1 + ' & resi ' + num_res1 + ' & n. CA'
    selection_string2 = 'chain' + chain_res2 + ' & resi ' + num_res2 + ' & n. CA'
    pymol.cmd.select(selection_string1 + ' ' + selection_string2)
    pymol.cmd.set_name('sele', pair_string)
    #label residues
    pymol.cmd.label(selection_string1,"'%s %s %s' %(resn,resi,chain)")
    pymol.cmd.label(selection_string2,"'%s %s %s' %(resn,resi,chain)")
    selected_residues = pairs_to_residues(selected_pairs)
    if pair[0] not in selected_residues:
        pymol.cmd.hide('labels', selection_string1)
    if pair[1] not in selected_residues:
        pymol.cmd.hide('labels', selection_string2)
    return
    
def pairs_to_residues(pairs):
    residues = []
    for i in range(len(pairs)):
        for j in range(len(pairs[0])):
            if pairs[i][j] not in residues:
                residues.append(pairs[i][j])
    return (residues)

def get_top_10(pairs, values):
    top_pairs = []
    absolute_values = []
    size_10_percent = len(values)//10
    for value in values:
        absolute_values.append(abs(value))
    absolute_values.sort(reverse=True)
    top_values = absolute_values[:size_10_percent]
    for f in range(len(pairs)):
        if len(top_pairs) <= len(top_values):
            if (values[f] in top_values) or (-1*values[f] in top_values):
                top_pairs.append(pairs[f])
    return (top_pairs)

def all_pairs_from_interest(pairs, residues_of_interest):
    selected_pairs = []
    for pair in pairs:
        if pair[0] in residues_of_interest or pair[1] in residues_of_interest:
            selected_pairs.append(pair)
    return (selected_pairs)

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
    
def generate_session(pdb_file, surfaces_file, residues_of_interest, color_scale_range, session_file_name):
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
    pairs, values = get_pairs_contacts(surfaces_file)
    if residues_of_interest is None:
        selected_pairs = get_top_10(pairs, values) # get pairs with largest absolute value of interaction - top 10%
    else:
        residues_of_interest = list(residues_of_interest[1:-1].split(","))
        selected_pairs = all_pairs_from_interest(pairs, residues_of_interest)
    color_codes = color_scale(values)
    for j in range(len(pairs)):
        color_distance(pairs[j], values[j], color_codes[j], selected_pairs)
        label_pairs(pairs[j], selected_pairs)
    show_separate_surfaces(chains)
    pymol.cmd.save(session_file_name, format='pse')
    return

### test function get_sum_per_residue
#residues, values = get_sum_per_residue('LIGIN_output_example.csv')
#print (color_scale(values))

### test function get_pairs_contacts
#pairs, values = get_pairs_contacts('output_7VQ0.csv')
#print (pairs)
#print (values)

### test function read_residue
#print (read_residue('TYR513B'))

def main():
    
    parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f","--pdb_file", action="store")
    parser.add_argument("-c","--input_csv_file", action="store")
    parser.add_argument("-o","--pymol_session_output_name", action="store")
    parser.add_argument("-cs","--color_scale_range", action="store")
    parser.add_argument("-res","--residues_of_interest", action="store")
    args=parser.parse_args()
    
    generate_session(args.pdb_file, args.input_csv_file, args.residues_of_interest, args.color_scale_range, args.pymol_session_output_name)
    
    return

main()
