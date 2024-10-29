# Load required libraries
library(graphite)
library(igraph)
library(org.Hs.eg.db)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide the KEGG BioServices result and an output directory.")
}

# Define input file and output directory from arguments
KEGG_Result_file <- args[1]
output_dir <- args[2]

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set up pathway and database
r <- pathways("hsapiens", "kegg")
hs <- org.Hs.eg.db

# Read the KEGG enrichment results
KEGG_Result <- read.csv(KEGG_Result_file, header = TRUE)
KEGG_Result$kegg_path <- gsub('hsa', 'hsa:', KEGG_Result$kegg_path)

# Process each pathway
for (i in 1:length(r)) {
  if (r[[i]]@id %in% KEGG_Result$kegg_path) {
    current.pathway <- convertIdentifiers(r[[i]], "SYMBOL")
    current.pathway_graphnel <- pathwayGraph(current.pathway, which = "mixed")
    current.pathway_graph <- graph_from_graphnel(current.pathway_graphnel, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
    
    # Clean up vertex names
    V(current.pathway_graph)$name <- gsub('SYMBOL:', '', V(current.pathway_graph)$name)
    
    # Convert to data frame and clean up columns
    current.pathway_data_frame <- as_long_data_frame(current.pathway_graph)
    current.pathway_data_frame <- current.pathway_data_frame[!apply(is.na(current.pathway_data_frame[, c(5, 6)]), 1, all), ]
    current.pathway_data_frame = as.data.frame(current.pathway_data_frame[c(5,6)])
    current.pathway_data_frame$ID <- gsub(' ', '_', r[[i]]@title)
    current_pathway_name <- gsub("[:/ ]", "_", paste0(r[[i]]@id, '_', r[[i]]@title, '_Kegg_edge_list.txt'))
    # Write to file in the specified output directory
    write.table(
      current.pathway_data_frame,
      file = file.path(output_dir, current_pathway_name),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }
}

