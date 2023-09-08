import pandas as pd

import argparse
import snakemake

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--outdir', metavar = '', type = str, default = "surfaces_output", 
    help = ""),
parser.add_argument('--snakefile', metavar = '', type = str,
    help = "")
parser.add_argument("-lig","--pdb_ligand", action="store", help="Expected in format PDBID_LIGANDCHAINID")
parser.add_argument("-dir","--pdb_dir", action="store")
parser.add_argument("-def","--atomtypes_definition", action="store")
parser.add_argument("-c","--chains", action="store")
parser.add_argument("-prefix","--def_prefix", action="store", default='custom', help = "prefix for the output def file")
parser.add_argument("-o","--output_name", action="store")
parser.add_argument("-dat","--atomtypes_interactions", action="store")
parser.add_argument("-vcon","--vcon_out", action="store")

args = parser.parse_args()



status = snakemake.snakemake(
    args.snakefile, 
    printshellcmds=False,
    configfiles={"config.yml"},
    cores = 1)