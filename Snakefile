import numpy as np
import pandas as pd
import os

from os import path

# Annotree beta uses GTDB 214.1
GTDB_VERSION = '214'
BASE_GTDB_URL = f"https://data.ace.uq.edu.au/public/gtdb/data/releases/release{GTDB_VERSION}/{GTDB_VERSION}.1/"
LOCAL_GTDB_DIR = "data/gtdb/"

# CD HIT parameters
CDHIT_SEQ_ID = 0.8
CDHIT_SUFFIX = f"{round(CDHIT_SEQ_ID,2)*100}"

# Current working directory
CWD = os.getcwd()
print('Current directory', CWD)

# We will fetch the following files from GTDB and unzip those that are gzipped
GTDB_FILENAMES2FETCH = [
    f"bac120_metadata_r{GTDB_VERSION}.tsv.gz",
    f"bac120_r{GTDB_VERSION}.tree",
    f"ar53_metadata_r{GTDB_VERSION}.tsv.gz",
    f"ar53_r{GTDB_VERSION}.tree",
]
GTDB_PATHS2FETCH = [BASE_GTDB_URL + n for n in GTDB_FILENAMES2FETCH]
LOCAL_GTDB_FNAMES = [path.join(LOCAL_GTDB_DIR, n) for n in GTDB_FILENAMES2FETCH]
LOCAL_GTDB_FNAMES = [n.strip('.gz') for n in LOCAL_GTDB_FNAMES]  # gonna unzip them
GTDB_BAC_METADATA, GTDB_BAC_TREE, GTDB_ARC_METADATA, GTDB_ARC_TREE = LOCAL_GTDB_FNAMES

annotree_manifest = "data/annotree/annotree_manifest.csv"
annotree_manifest_df = pd.read_csv(annotree_manifest).dropna(how='all')
ANNOTREE_QS = tuple(annotree_manifest_df["query"].unique())
ANNOTREE_NAMES = tuple(annotree_manifest_df["name"].unique())
print(annotree_manifest_df)
print(ANNOTREE_QS)
ANNOTREE_Q_MAP = annotree_manifest_df.set_index('query')['name'].to_dict()
ANNOTREE_Q_NAMES = [ANNOTREE_Q_MAP[q].replace(' ', '_') for q in ANNOTREE_QS]
ANNOTREE_Q_NAMES_STR = " ".join(ANNOTREE_Q_NAMES)

annotree_manifest_df = annotree_manifest_df.set_index("domain,query".split(','))

NUTRIENTS = annotree_manifest_df.nutrient.unique()


# Demands final output is compiled
rule all:
    input:
        gtdb_data=LOCAL_GTDB_FNAMES,
        gtdb_stats="output/gtdb_phylo_stats.csv",
        pi_funcs_by_organism="output/gene_funcs_by_organism.csv",
        itol_bac_tree=expand("output/itol_bac_{nutrient}_phylum.txt", nutrient=NUTRIENTS)

# Fetch the GTDB metadata and tree for the version specified
rule fetch_gtdb:
    output:
        LOCAL_GTDB_FNAMES
    shell:
        "mkdir -p {LOCAL_GTDB_DIR} && "
        "cd {LOCAL_GTDB_DIR} && "
        # need -O for each file to avoid overwriting
        "for url in {GTDB_PATHS2FETCH}; do curl -O \"$url\"; done && "
        "gunzip -f *.gz"

rule calc_gtdb_stats:
    input:
        bacteria=GTDB_BAC_METADATA,
        archaea=GTDB_ARC_METADATA
    output:
        "output/gtdb_phylo_stats.csv"
    shell:
        "python scripts/gtdb2stats.py --representatives_only -b {input.bacteria} -a {input.archaea} -o {output}"

rule calc_pi_funcs_by_organism:
    input:
        manifest=annotree_manifest,
    output:
        outdir="intermediate/annotree/",
        wide_output="genes_by_organism.csv",
    shell:
        "python scripts/tabulate_gene_funcs.py --manifest {input.manifest} "
        "--outdir {output.outdir} --outfname {output.wide_output}"

rule apply_boolean_expressions:
    input:
        genes_by_organism="intermediate/annotree/genes_by_organism.csv",
        expressions_fname="data/annotree/annotree_expressions.csv"
    output:
        outdir="intermediate/annotree/",
        nutrient_outputs=expand("intermediate/annotree/{nutrient}_functional_results.csv", nutrient=NUTRIENTS)
    shell:
        "python scripts/apply_expressions.py --input {input.genes_by_organism} --expressions {input.expressions_fname} --outdir {output.outdir}"

rule make_itol_tree:
    input: 
        "intermediate/annotree/{nutrient}_functional_results.csv"
    output:
        "output/itol_bac_{nutrient}_phylum.txt"
    shell:
        "python scripts/hits2itol.py --in {input} --out {output} -d bacteria --agg_level phylum -t heatmap"

