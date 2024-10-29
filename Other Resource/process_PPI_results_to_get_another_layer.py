import pandas as pd
import argparse

# Define a function to fetch unique genes from a given input file
def get_PPI_genes(input_file):
    df = pd.read_csv(input_file, sep='\t')
    df1 = df[['Official Symbol Interactor A', 'Official Symbol Interactor B']]
    unique_genes = pd.DataFrame(set(df1['Official Symbol Interactor A']) | set(df1['Official Symbol Interactor B']))
    unique_genes = "'" + unique_genes + "'"
    unique_genes_T = unique_genes.T
    unique_genes_T['all'] = unique_genes_T.apply(lambda x: ', '.join(x[x.notnull()]), axis=1)
    unique_genes_T = unique_genes_T.astype(str)
    unique_genes_T_1 = unique_genes_T['all']
    return unique_genes_T_1

# Parse command-line arguments for input and output files
parser = argparse.ArgumentParser(description="Process gene interactions for a given input file and save unique genes to an output file.")
parser.add_argument("input_file", help="Path to the input file containing gene interactions")
parser.add_argument("output_file", help="Path to the output file where unique genes will be saved")
args = parser.parse_args()

# Process the input file and save results to the output file
output = get_PPI_genes(args.input_file)
output.to_csv(args.output_file, index=None, header=None, sep='\t')

