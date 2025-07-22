import pandas as pd
import argparse
import os

from os import path


def main():
    parser = argparse.ArgumentParser(description="Calculate PI functions by organism from the manifest file.")
    parser.add_argument('--manifest', type=str, required=True, help='Path to the annotree manifest file')
    parser.add_argument('--out_long', type=str, default='genes_long.csv',
                        help='Path to the long-format output file.')
    parser.add_argument('--out_wide', type=str, default='genes_by_organism.csv',
                        help='Path to the wide-format output file.')
    args = parser.parse_args()

    outdir = path.dirname(args.out_wide)
    print(f"Will write ancillary files to {outdir}")
    if not path.exists(outdir):
        print(f"Creating output directory: {outdir}")
        os.makedirs(outdir, exist_ok=True)

    manifest = pd.read_csv(args.manifest).dropna(how='all')
    manifest['filename'] = manifest['name'].apply(lambda x: x + '_bacteria.csv')
    manifest['filepath'] = manifest['filename'].apply(lambda x: f"data/annotree/{x}")

    # simply combine all the dataframes into a single one with a common index
    all_results = []
    for _, row in manifest.iterrows():
        df = pd.read_csv(row['filepath'], index_col=0)
        all_results.append(df)

    # concatenate all results into a single DataFrame along the rows
    combined_df = pd.concat(all_results, axis=0)
    combined_df = combined_df.reset_index()

    # check for duplicate geneId values
    duplicates = combined_df.duplicated(subset=['geneId'], keep=False)
    if duplicates.any():
        print(f"found {duplicates.sum()} duplicate gene IDs in a total of {len(combined_df)} genes.")
        dups = combined_df[duplicates]
        dup_queries = dups['SearchId'].unique()
        print(f"Duplicates from queries: {dup_queries}")
        dup_q_functions = manifest.set_index('name').loc[dup_queries]['function']
        print(f"Duplicate query gene functions: {dup_q_functions}")
        
        # Save the duplicates for further inspection
        duplicates_file = path.join(outdir, 'gene_duplicates_long.csv')
        dups.to_csv(duplicates_file, index=False)
        print(f"Duplicates saved to {duplicates_file}")

        # We are keeping the duplicates since they may represent gene 
        # families that are homologous but perform different functions.
        # Hard to know which one to keep, so we keep them all.

    # Save the combined DataFrame to the specified output directory
    combined_df.to_csv(args.out_long, index=False)
    print(f"Combined data saved to {args.out_long}")

    # Make a wide format DataFrame with gtdbId (species) as index and SearchId as columns
    # values are boolean indicating presence of the gene

    # Drop duplicates first to ensure unique index/column pairs for pivot
    # This is fine because we are binarizing the presence of genes in the end anyway
    combined_df_nodup = combined_df.drop_duplicates(subset=['gtdbId', 'SearchId'])
    wide_df = combined_df_nodup.pivot(index='gtdbId', columns='SearchId', values='geneId').notnull()
    wide_df.reset_index(inplace=True)

    # Save the wide format DataFrame to the specified output directory
    wide_df.to_csv(args.out_wide, index=False)
    print(f"Wide format data saved to {args.out_wide}")


if __name__ == "__main__":
    main()