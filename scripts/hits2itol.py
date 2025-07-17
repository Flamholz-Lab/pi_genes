#!/usr/bin/env python

"""Script to convert AnnoTree hits to iTOL dataset.

TODO: automate legend generation.
TODO: make sure that the AnnoTree tree is using the GTDB taxonomy.
"""

__author__ = 'Avi I. Flamholz'

import argparse
import pandas as pd
import re
import seaborn as sns

from os import path
from pathlib import Path

REPS_FNAMES = {
    'bacteria': 'bac120_metadata_r214.tsv',
    'archaea': 'ar53_metadata_r214.tsv'
}

SELF_PATH = Path(path.abspath(__file__))
GTDB_PATH = path.join(SELF_PATH.parent.parent, 'data/gtdb/')

ANNOTREE_ID = "gtdb_id"
SIMPLEBAR_HEADER_FORMAT = """
DATASET_SIMPLEBAR
SEPARATOR COMMA
DATASET_LABEL,{label}
COLOR,{hex_color}
DATA
"""

MULTIBAR_HEADER_FORMAT = """
DATASET_MULTIBAR
SEPARATOR COMMA
DATASET_LABEL,{label}
COLOR,#ff0000
FIELD_COLORS,{field_colors}
FIELD_LABELS,{field_labels}
DATA
"""

HEATMAP_HEADER_FORMAT = """
DATASET_HEATMAP
SEPARATOR COMMA
DATASET_LABEL,{label}
COLOR,#000000
COLOR_NAN,#000000
LEGEND_COLORS,{field_colors}
LEGEND_LABELS,{field_labels}
LEGEND_SHAPES,{field_shapes}
FIELD_COLORS,{field_colors}
FIELD_LABELS,{field_labels}
FIELD_SHAPES,{field_shapes}
DATA
"""

BINARY_HEADER_FORMAT = """
DATASET_BINARY
SEPARATOR COMMA
DATASET_LABEL,{label}
COLOR,#ff0000
FIELD_COLORS,{field_colors}
FIELD_LABELS,{field_labels}
FIELD_SHAPES,{field_shapes}
DATA
"""

# A regex capturing all the taxonomic levels in the GTDB taxonomy format
TAXONOMY_RE = re.compile(r'd__(?P<domain>[^;]+);p__(?P<phylum>[^;]+);c__(?P<class>[^;]+);o__(?P<order>[^;]+);f__(?P<family>[^;]+);g__(?P<genus>[^;]+);s__(?P<species>[^;]+)')

AGG_LEVELS = {
    'species': 'gtdbId',
    'genus': 'genus',
    'family': 'family',
    'order': 'order',
    'class': 'class',
    'phylum': 'phylum'
}
 
def count_hits(gids, reps_df, agg_level):
    """Count hits in annotree CSV file at this aggregation level.
    
    Args:
        fh (file): File handle to annotree CSV file.
        agg_level (str): Aggregation level to count hits at.
    
    Returns:
        pd.DataFrame: DataFrame with counts at this aggregation level.
    """
    masked_reps_df = reps_df[reps_df.index.isin(gids)]
    agg_key = AGG_LEVELS[agg_level]
    counts = masked_reps_df.groupby(agg_key).agg(dict(species='count'))
    counts.columns = ['count']
    return counts


def normalize_counts(counts, reps_df, agg_level):
    """Normalize hit counts by number of representative genomes.
    
    Args:
        counts (pd.DataFrame): DataFrame with counts at this aggregation level.
        reps_df (pd.DataFrame): DataFrame with representative genome taxonomy.
        agg_level (str): Aggregation level to normalize counts at.
    
    Returns:
        pd.DataFrame: DataFrame with normalized counts at this aggregation level.
    """
    # No normalization needed for species level
    if agg_level == 'species':
        return counts
    
    # Count number of representative genomes at this aggregation level
    reps_count = reps_df.groupby(agg_level).agg(dict(species='count'))
    reps_count.columns = ['count']

    # Normalize counts by number of representative genomes
    normed = counts / reps_count
    mask = normed['count'].notnull()

    # Drop rows with no data.
    return normed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--domain', type=str, default='bacteria',
                        help='Domain to use for GTDB representative genomes.',
                        choices=('bacteria', 'archaea'))
    parser.add_argument('--in', '-i', type=argparse.FileType('r'),
                        required=True, dest='input',
                        help='Input annotree hits CSV file path. Indexed by gtdbId.')
    parser.add_argument('--out', '-o', type=str, default='iTOL_dataset.txt',
                        help='Output iTOL dataset file name.')
    parser.add_argument('--palette', '-p', type=str, default='tab10',
                        help='Seaborn palette name for colors.')
    parser.add_argument('--agg_level', '-a', type=str, default='species',
                        help='Aggregation level for annotree hits.',
                        choices=sorted(AGG_LEVELS.values()))
    parser.add_argument('--plot_type', '-t', type=str, default='bar',
                        help='iTOL plot type.',
                        choices=('bar', 'binary', 'heatmap'))
    parser.add_argument('--binary_threshold', '-b', type=float, default=0.5,
                        help='Threshold for binarizing counts or normalized counts.')
    
    args = parser.parse_args()
    print(f'Input file: {args.input.name}')

    print('Reading representatives...')
    gtdb_reps_fname = path.join(GTDB_PATH, REPS_FNAMES[args.domain])
    reps_df = pd.read_csv(gtdb_reps_fname, sep='\t')

    # Extract the taxonomic information 
    tax_df = reps_df['gtdb_taxonomy'].str.extract(TAXONOMY_RE)
    reps_df = pd.concat([reps_df, tax_df], axis=1)

    # Retain only representative genomes
    mask = reps_df['gtdb_representative'] == 't'
    reps_df = reps_df[mask].set_index('accession')

    hits = pd.read_csv(args.input, index_col=0)
    count_cols = []
    for c in hits.columns:
        print('Processing column:', c)
        gids = hits[hits[c] == True].index.to_list()
        counts = count_hits(gids, reps_df, args.agg_level)
        normed = normalize_counts(counts, reps_df, args.agg_level)
        normed.columns = [c]
        count_cols.append(normed)
    annotree_counts = pd.concat(count_cols, axis=1)

    print('Aggregating counts...')
    labels = annotree_counts.columns.tolist()
    n_colors = len(labels)
    pal = sns.color_palette(args.palette, n_colors=n_colors)
    hex_colors = pal.as_hex()

    header_text = ''
    if args.plot_type == 'binary':
        print('Binarizing counts...')
        field_colors = ','.join(hex_colors)
        field_labels = ','.join(labels)
        field_shapes = ','.join(['1' for _ in range(n_colors)])
        header_text = BINARY_HEADER_FORMAT.format(
            label='annotree', field_colors=field_colors,
            field_labels=field_labels, field_shapes=field_shapes)
        annotree_counts = (annotree_counts > args.binary_threshold).astype(int)
    elif args.plot_type == 'bar' and n_colors > 1:
        print('Creating multibar plot...')
        field_colors = ','.join(hex_colors)
        field_labels = ','.join(labels)
        header_text = MULTIBAR_HEADER_FORMAT.format(
            label='annotree', field_colors=field_colors,
            field_labels=field_labels)
    elif args.plot_type == 'heatmap':
        print('Creating heatmap plot...')
        field_colors = ','.join(hex_colors)
        field_labels = ','.join(labels)
        field_shapes = ','.join(['1' for _ in range(n_colors)])
        header_text = HEATMAP_HEADER_FORMAT.format(
            label='annotree', field_colors=field_colors,
            field_labels=field_labels, field_shapes=field_shapes)
    else:
        print('Creating simplebar plot...')
        header_text = SIMPLEBAR_HEADER_FORMAT.format(
            label=args.labels[0], hex_color=hex_colors[0])
    
    print(f'Writing iTOL dataset to {args.out}...')
    with open(args.out, 'w') as f:
        f.write(header_text)
        annotree_counts.to_csv(f, header=False)
    
    print('Done!')


if __name__ == '__main__':
    main()
