from exp_parsing import BooleanExpressionParser
import pandas as pd
import argparse
from os import path


def main():
    parser = argparse.ArgumentParser(description="Apply boolean expressions to gene rows.")
    parser.add_argument('--input', type=str, required=True, help='Path to the input CSV file with gene data. Rows represent genomes and columns represent genes.')
    parser.add_argument('--expressions', type=str, required=True, help='Path to the file containing boolean expressions')
    parser.add_argument('--outdir', type=str, required=True, help='Path to the output directory')
    args = parser.parse_args()

    # Load the gene data
    gene_data_df = pd.read_csv(args.input, index_col=0).dropna(how='all')

    # Load the expressions -- rows are functional categories with expressions to be 
    # run in the same order as the file
    expressions_df = pd.read_csv(args.expressions, index_col=0).dropna(how='all')

    parser = BooleanExpressionParser()
    for function_name, row in expressions_df.iterrows():
        parsed = parser.parse_expression(row['boolean_expression'])
        print(f"Applying expression for {function_name}: {row['boolean_expression']}")
        print(f"Parsed expression: {parsed.as_list()}")

        gene_data_df[function_name] = gene_data_df.apply(
            lambda row: parser.evaluate(parsed, row), axis=1
        )

    # Save the results to the output file
    full_output_path = path.join(args.outdir, 'gene_data_with_derived_functions.csv')
    gene_data_df.to_csv(full_output_path)
    print(f"Results saved to {full_output_path}")

    # Retain only the columns marked as "for_display" in the expressions DataFrame
    display_columns = expressions_df[expressions_df['for_display']].index.tolist()
    display_df = gene_data_df[display_columns]
    display_output_path = path.join(args.outdir, 'functional_results.csv')
    display_df.to_csv(display_output_path)
    print(f"Display results saved to {display_output_path}")

    expressions_for_display = expressions_df[expressions_df['for_display']]
    for nutrient in expressions_for_display['nutrient'].unique():
        nutrient_exps = expressions_for_display[expressions_for_display['nutrient'] == nutrient].index.tolist()
        nutrient_df = display_df[nutrient_exps]

        nutrient_output_path = path.join(args.outdir, f'{nutrient}_functional_results.csv')
        nutrient_df.to_csv(nutrient_output_path)
        print(f"Results for {nutrient}-related genetic functions saved to {nutrient_output_path}")

if __name__ == "__main__":
    main()