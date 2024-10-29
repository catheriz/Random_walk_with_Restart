from bioservices import KEGG
import pandas as pd
import argparse

def fetch_kegg_pathways(disease_gene_file, output_file):
    k = KEGG(verbose=True)
    k.organism = "hsa"
    s = KEGG()

    # Read disease genes from the input file
    with open(disease_gene_file, 'r') as file:
        disease_genes = [line.strip() for line in file if line.strip()]
    print("Disease genes loaded from file:", disease_genes)
    # Initialize lists to store results
    Kegg_path = []
    Kegg_name = []
    gene_name = []

    # Fetch KEGG pathways for each gene
    for gene in disease_genes:
        try:
            for path, name in s.get_pathway_by_gene(gene, "hsa").items():
                gene_name.append(gene)
                Kegg_path.append(path)
                Kegg_name.append(name)
        except (AttributeError):
            pass

    # Create DataFrame and save to CSV
    KEGG_pathway = pd.DataFrame({'kegg_name': Kegg_name, 'kegg_path': Kegg_path, 'gene_name': gene_name})
    KEGG_pathway.to_csv(output_file, header=True, index=False)

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Fetch KEGG pathways for a list of disease genes.")
parser.add_argument("disease_gene_file", help="Path to the input file containing a list of disease genes.")
parser.add_argument("output_file", help="Output file name for the KEGG pathway data.")
args = parser.parse_args()

# Run the function
fetch_kegg_pathways(args.disease_gene_file, args.output_file)

