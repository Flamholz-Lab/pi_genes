import pandas as pd
import argparse

from os import path


def main():
    parser = argparse.ArgumentParser(description="Calculate PI functions by organism from the manifest file.")
    parser.add_argument('--manifest', type=str, required=True, help='Path to the annotree manifest file')
    parser.add_argument('--output', type=str, required=True)
    parser.add_argument('--outdir', type=str, default='intermediate/annotree/',
                        help='Output directory for the combined data')
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest)
    manifest['filename'] = manifest['name'].apply(lambda x: x + '_bacteria.csv')
    manifest['filepath'] = manifest['filename'].apply(lambda x: f"data/annotree/{x}")

    all_results_by_ft = dict()
    for nutrient, nutrient_df in manifest.groupby('nutrient'):
        results_dfs_by_ft = dict()
        print(f"Processing gene families associated with {nutrient}")
        for gname, group_df in nutrient_df.groupby('function_type'):
            row_dfs = []
            for idx, row in group_df.iterrows():
                print(f"Processing {row['name']} with function type {row['function_type']}")
                row_dfs.append(pd.read_csv(row['filepath']))
            combined_df = pd.concat(row_dfs, ignore_index=True)
            
            output_path = path.join(args.outdir, f"{gname}_bacteria.csv")
            combined_df.to_csv(output_path, index=False)
            print(f"Saved combined data for {gname} to {output_path}")

            unique_organisms = combined_df['gtdbId'].unique()
            print(f"Unique organisms for {gname}: {len(unique_organisms)}")
            results_dfs_by_ft[gname] = unique_organisms
            all_results_by_ft[gname] = unique_organisms

        all_organism_ids = set()
        for org_ids in results_dfs_by_ft.values():
            all_organism_ids.update(org_ids)
        all_organism_ids = list(all_organism_ids)  # Convert to list for DataFrame creation
        print(f"Total unique organisms across all function types: {len(all_organism_ids)}")

        # Create a DataFrame with all unique organisms
        func_df = pd.DataFrame(index=all_organism_ids)
        func_df.index.name = 'gtdbId'
        for ft in results_dfs_by_ft.keys():
            func_df[ft] = func_df.index.isin(results_dfs_by_ft[ft])

        out_fname = f"{nutrient}_gene_funcs_by_organism.csv"
        out_path = path.join(args.outdir, out_fname)
        func_df.to_csv(out_path, index=True)
        print(f"Saved combined PI functions by organism to {out_fname}")
        print(func_df.head())

    # Combine all function types into a single DataFrame
    all_funcs_df = pd.DataFrame(index=all_results_by_ft.keys())
    all_funcs_df.index.name = 'gtdbId'
    for ft, org_ids in all_results_by_ft.items():
        all_funcs_df[ft] = all_funcs_df.index.isin(org_ids)
    all_funcs_df.to_csv(args.output, index=True)
    print(f"Saved all PI functions by organism to {args.output}")

if __name__ == "__main__":
    main()