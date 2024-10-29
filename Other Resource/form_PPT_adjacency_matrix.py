import pandas as pd
import numpy as np
import networkx as nx
import argparse

def make_PPI_matrix(layer_files, combined_output, adjacency_output):
    all_data = []  # Collect all layers’ data
    unique_genes_sets = []  # Collect unique genes from each layer

    # Process each layer file provided by the user
    for idx, layer_file in enumerate(layer_files):
        df_layer = pd.read_csv(layer_file, sep='\t')  # Read each layer file
        if 'Official Symbol Interactor A' not in df_layer.columns or 'Official Symbol Interactor B' not in df_layer.columns:
            raise ValueError(f"The file {layer_file} is missing required columns 'Official Symbol Interactor A' or 'Official Symbol Interactor B'.")

        # Extract unique genes in this layer
        unique_genes = list(set(df_layer['Official Symbol Interactor A']) | set(df_layer['Official Symbol Interactor B']))
        
        # Save unique genes for each layer to individual output files
        unique_genes_df = pd.DataFrame(unique_genes, columns=['Gene'])
        unique_genes_df.to_csv(f"Layer{idx+1}_unique_genes.csv", header=True, index=False)
        
        # Collect data for concatenation and unique genes
        all_data.append(df_layer)
        unique_genes_sets.append(set(unique_genes))

    # Combine data from all layers into a single DataFrame
    combined_data = pd.concat(all_data, axis=0)
    combined_genes = pd.DataFrame(
        np.sort(combined_data[['Official Symbol Interactor A', 'Official Symbol Interactor B']]), 
        columns=['Official Symbol Interactor A', 'Official Symbol Interactor B']
    ).value_counts().reset_index(name='counts')
    
    # Save combined edge list with weights
    combined_genes.to_csv(combined_output, header=True, index=True)
    
    # Generate adjacency matrix
    G = nx.from_pandas_edgelist(combined_genes, 'Official Symbol Interactor A', 'Official Symbol Interactor B', edge_attr='counts')
    adj = nx.adjacency_matrix(G, weight='counts')
    adj_matrix = pd.DataFrame(adj.todense(), index=G.nodes, columns=G.nodes)
    
    # Save the adjacency matrix
    adj_matrix.to_csv(adjacency_output, header=True, index=True)

    return adj_matrix

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process one or more gene interaction layers and generate a PPI adjacency matrix.")
parser.add_argument("layer_files", nargs='+', help="Paths to the gene interaction files for each layer (at least one required)")
parser.add_argument("combined_output", help="Output path for the combined gene edge list")
parser.add_argument("adjacency_output", help="Output path for the adjacency matrix CSV file")
args = parser.parse_args()

# Ensure at least one layer file is provided
if len(args.layer_files) < 1:
    raise ValueError("Please provide at least one layer file.")

# Generate and save the adjacency matrix
adj_matrix = make_PPI_matrix(args.layer_files, args.combined_output, args.adjacency_output)

