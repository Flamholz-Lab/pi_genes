# a script that takes two gtdb metadata files as input,
# one for bacteria and one for archaea. the script then 
# counts up the number of genomes at each phylogenetic level

import pandas as pd
import argparse

from os import path

# Full names of the phylogenetic levels
FULL_PHYLO_NAMES_DICT = {
    'd': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'
}
PHYLO_COLNAMES = 'domain,phylum,class,order,family,genus,species'.split(',')

# Name of the column in the annotree results that contains the concatenated phylogeny
TAXONOMY_COLNAME = 'gtdb_taxonomy'

def extract_phylogeny(row):
    # phylogeny is given as a string with the following example
    # d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosopumilus;s__Nitrosopumilus sp905612295
    # extract a dictionary mapping levels to values
    phylo = {}
    for part in row[TAXONOMY_COLNAME].split(';'):
        key, value = part.split('__')
        # Use the full phylogeny names
        phylo[FULL_PHYLO_NAMES_DICT[key]] = value
    return phylo


def add_phylogeny_columns(df):
    # apply function to each row in the dataframe
    phylo_df = df.apply(extract_phylogeny, axis=1, result_type='expand')
    return pd.concat([df, phylo_df], axis=1)

def metadata_phylogeny_counts(metadata_df):
    all_counts = []

    my_df = add_phylogeny_columns(metadata_df)

    # for each column in PHYLO_COLNAMES, we want counts of values by name
    cts_by_phylo = dict(phylogenetic_level=[], domain=[], name=[], count=[])
    for colname in PHYLO_COLNAMES:
        my_cts = my_df.groupby('domain')[colname].value_counts().reset_index()
        for _, row in my_cts.iterrows():
            cts_by_phylo['phylogenetic_level'].append(colname)
            cts_by_phylo['domain'].append(row['domain'])
            cts_by_phylo['name'].append(row[colname])
            cts_by_phylo['count'].append(row['count'])

    # convert to a dataframe
    cts_by_phylo_df = pd.DataFrame(cts_by_phylo)

    # convert to long format, rename index to 'name' and variable to 'phylogenetic_level'
    #cts_by_phylo_df = pd.melt(cts_by_phylo_df, id_vars='index', value_vars=PHYLO_COLNAMES,
    #                            var_name='phylogenetic_level', value_name='count')
    #cts_by_phylo_df = cts_by_phylo_df.rename(columns={'index': 'name'})

    # Drop rows with NaN values
    #cts_by_phylo_df = cts_by_phylo_df.dropna()
    return cts_by_phylo_df

def main():
    parser = argparse.ArgumentParser(description='Calculate phylogenetic summary statistics from a manifest file')
    parser.add_argument('-b', '--bacterial_metadata', type=str, help='Input bacterial metadata file')
    parser.add_argument('-a', '--archaeal_metadata', type=str, help='Input archaeal metadata file')
    parser.add_argument('-o', '--output', type=str, help='Output summary statistics file name')
    parser.add_argument('-r', '--representatives_only', action='store_true',
                        help='Limit to representative genomes only',
                        default=False)
    parser.add_argument('--sep', type=str, help='Default separator is a tab', default='\t')

    args = parser.parse_args()

    print(f'Reading bacterial metadata from {args.bacterial_metadata}')
    print(f'Reading archaeal metadata from {args.archaeal_metadata}')
    bacterial_df = pd.read_csv(args.bacterial_metadata, sep=args.sep)
    archaeal_df = pd.read_csv(args.archaeal_metadata, sep=args.sep)

    if args.representatives_only:
        bacterial_df = bacterial_df[bacterial_df['gtdb_representative'] == 't']
        archaeal_df = archaeal_df[archaeal_df['gtdb_representative'] == 't']

    bacterial_counts = metadata_phylogeny_counts(bacterial_df)
    archaeal_counts = metadata_phylogeny_counts(archaeal_df)

    all_counts_df = pd.concat([bacterial_counts, archaeal_counts], axis=0)

    print(f'Writing summary statistics to {args.output}')
    all_counts_df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()