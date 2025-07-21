import pandas as pd
import argparse

from os import path


def main():
    parser = argparse.ArgumentParser(description="Calculate PI functions by organism from the manifest file.")
    parser.add_argument('--manifest', type=str, required=True, help='Path to the annotree manifest file')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Path to the output directory')
    parser.add_argument('--outfname', type=str, default='genes_by_organism.csv',
                        help='Name of the wide-format output file')
    args = parser.parse_args()

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
        
        # Save the file with duplicates for further inspection
        duplicates_file = path.join(args.outdir, 'results_long_duplicates.csv')
        dups.to_csv(duplicates_file, index=False)
        print(f"Long-form results with duplicates saved to {duplicates_file}")

        # Drop the duplicates -- TODO: do smarter here.
        print(f"Dropping {duplicates.sum()} duplicates, keeping the first occurrence.")
        combined_df = combined_df.drop_duplicates(subset=['geneId'], keep='first')
        combined_df = combined_df.reset_index()

    # Save the combined DataFrame to the specified output directory
    output_file_long = path.join(args.outdir, 'results_combined.csv')
    combined_df.to_csv(output_file_long, index=False)
    print(f"Combined data saved to {output_file_long}")

    # Make a wide format DataFrame with gtdbId (species) as index and SearchId as columns
    # values are boolean indicating presence of the gene

    # Drop duplicates first to ensure unique index/column pairs for pivot
    # This is fine because we are binarizing the presence of genes in the end anyway
    combined_df_nodup = combined_df.drop_duplicates(subset=['gtdbId', 'SearchId'])
    wide_df = combined_df_nodup.pivot(index='gtdbId', columns='SearchId', values='geneId').notnull()
    wide_df.reset_index(inplace=True)

    # Save the wide format DataFrame to the specified output directory
    output_file_wide = path.join(args.outdir, args.outfname)
    wide_df.to_csv(output_file_wide, index=False)
    print(f"Wide format data saved to {output_file_wide}")


if __name__ == "__main__":
    main()