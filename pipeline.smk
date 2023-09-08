pepfile: "/users/matthewcrown/GitHub/surfaces/config.yml"
PDB_LIGANDS = pep.sample_table["sample_name"].unique()
CIF_IDS = pep.sample_table["pdb_id"].unique()
print(PDB_LIGANDS)
rule all: 
    input:
        config['outdir'] + "/all_pdb_contacts.pkl",
        expand(config['outdir'] + "/pdb/{pdb_id}.pdb", pdb_id = CIF_IDS)        
rule collate_contacts:
    input:
        expand(config['outdir'] + "/surfaces/List_{pdb_ligand}.txt", pdb_ligand = PDB_LIGANDS)
    output:
        config['outdir'] + "/all_pdb_contacts.pkl"
    shell:
        "touch {config[outdir]}/all_pdb_contacts.pkl"#"collate_contacts.py -i {config[outdir]}/contacts -o {output}"

rule run_surfaces:
    input:
        pdb_file = config['indir'] + "/surfaces_cath_residue_df_copy.csv"
    output:
        matrix_file = config['outdir'] + "/surfaces/{pdb_ligand}.csv",
        list_file = config['outdir'] + "/surfaces/List_{pdb_ligand}.txt"
    params:
        pdb_dir = config['outdir'] + "/pdb/"
    shell:
        "python3 ligands_surfaces_all.py -lig_id {wildcards.pdb_ligand} -lig_file {input} -dir {params.pdb_dir} -def AMINO_FlexAID.def -dat FlexAID.dat -prefix {wildcards.pdb_ligand}_prefix -o {wildcards.pdb_ligand} -outdir {config[outdir]}/surfaces -vcon {wildcards.pdb_ligand}_Vcon_out.txt"

rule ciftopdb:
    input:
        cif_file = config['indir'] + "/cif/{pdb_id}.cif"
    output:
        pdb_file = config['outdir'] + "/pdb/{pdb_id}.pdb"
    shell:
        "python3 cif2pdb.py {input.cif_file} {output.pdb_file}"