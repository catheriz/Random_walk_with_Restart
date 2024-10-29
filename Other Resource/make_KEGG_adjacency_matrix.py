import pandas as pd
import networkx as nx
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process KEGG interaction data and output adjacency matrix.")
parser.add_argument("input_file", help="Path to the input file containing KEGG edge list data.")
parser.add_argument("output_edge_list", help="Path for saving the processed edge list with weights.")
parser.add_argument("output_adj_matrix", help="Path for saving the adjacency matrix.")
args = parser.parse_args()

# Load the KEGG edge list data
Kegg = pd.read_csv(args.input_file, sep=' ', header=None, names=['colA', 'colB', 'ID'])

# Filter out rows where 'colA' or 'colB' contain ":"
Kegg = Kegg[~(Kegg.colA.str.contains(":") | Kegg.colB.str.contains(":"))].reset_index(drop=True)

# Process directed edge list with weights
Kegg_genes = pd.DataFrame(Kegg.groupby(['colA', 'colB']).size()).reset_index()
Kegg_genes.rename(columns={0: 'counts'}, inplace=True)
Kegg_genes['Type'] = 'KEGG'

# Save the directed edge list with weights
Kegg_genes.to_csv(args.output_edge_list, header=True, index=False)

# Create directed graph and remove self-loops
Kegg_G = nx.from_pandas_edgelist(Kegg_genes, source='colA', target='colB', create_using=nx.DiGraph, edge_attr=['counts', 'Type'])
Kegg_G.remove_edges_from(nx.selfloop_edges(Kegg_G))

# Generate and save adjacency matrix
Kegg_adj = nx.to_pandas_adjacency(Kegg_G, weight='counts')
Kegg_adj.to_csv(args.output_adj_matrix, header=True, index=True)

