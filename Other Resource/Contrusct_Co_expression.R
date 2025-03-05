# Load necessary libraries
library(reshape2)
library(igraph)
library(biomaRt)

# Read tissue and immune cell expression data
data.tissue <- read.csv('RNA_consensus_tissue_gene_data/rna_tissue_consensus.tsv', head = TRUE, sep = '\t')
data.celline <- read.csv('RNA_HPA_immune_cell_sample_gene_data/rna_immune_cell_sample.tsv', head = TRUE, sep = '\t')

# Select relevant columns
data.celline <- data.celline[c(4,5,3,8)]

# Standardize column names
colnames(data.tissue) <- c('Ensemble_ID', 'Gene.name', 'Sample', 'Value')
colnames(data.celline) <- c('Ensemble_ID', 'Gene.name', 'Sample', 'Value')
# Combine datasets
data_Coexpr <- rbind(data.tissue, data.celline)

# Function to check and validate HGNC symbols
check.HGNC.symbols <- function(GeneNames) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensemble_HGNC <- getBM(attributes = c('hgnc_symbol'), 
                         filters = 'hgnc_symbol', 
                         values = GeneNames, 
                         mart = ensembl)
  return(ensemble_HGNC$hgnc_symbol)
}

# Validate HGNC gene names
ensemble_HGNC <- check.HGNC.symbols(unique(data_Coexpr$Gene.name))
data3_Coexpr <- data_Coexpr[data_Coexpr$Gene.name %in% ensemble_HGNC, ]

# Create expression matrix
matrix_Coexpr <- dcast(data3_Coexpr, Gene.name ~ Sample, mean, value.var = "Value")

# Set row names
rownames(matrix_Coexpr) <- matrix_Coexpr[, 1]
matrix_Coexpr <- matrix_Coexpr[, -1]

# Remove rows with missing values
matrix_Coexpr <- na.omit(matrix_Coexpr)

# Identify rows and columns with zero standard deviation
zero_sd_rows <- apply(matrix_Coexpr, 1, sd) == 0
zero_sd_cols <- apply(matrix_Coexpr, 2, sd) == 0

# Remove zero-variance rows and columns
matrix_Coexpr <- matrix_Coexpr[!zero_sd_rows, ]
matrix_Coexpr <- matrix_Coexpr[, !zero_sd_cols]

# Compute Spearman correlation matrix
CR <- cor(t(matrix_Coexpr), method = "spearman")

# Set diagonal to zero (self-correlations are not meaningful)
diag(CR) <- 0

# Extract significant correlations (absolute value ≥ 0.75)
inds <- which(abs(CR) >= 0.75, arr.ind = TRUE)
inds <- as.data.frame(inds)

# Extract correlation values
inds$weight <- abs(CR[inds$row, inds$col])

# Retrieve gene names
inds$row_name <- rownames(CR)[inds$row]
inds$col_name <- colnames(CR)[inds$col]

# Normalize network weights
inds$network_weight <- inds$weight / min(inds$weight)

# Create edge list for network analysis
co_express_network_edge_list <- inds[, c("row_name", "col_name", "network_weight")]

# Save edge list to CSV
write.csv(co_express_network_edge_list, 'Coexpression_network_edgelist_corr0.75_with_weight.csv', 
          quote = FALSE, row.names = FALSE)
