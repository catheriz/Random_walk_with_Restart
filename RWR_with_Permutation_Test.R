#!/usr/bin/env Rscript

# Load necessary libraries
library(igraph)
library(ggplot2)
library(dplyr)
library(Matrix)
library(Rcpp)
library(biomaRt)
library(supraHex)
library(dnet)
library(Rgraphviz)
library(data.table)

# Set seed for reproducibility
set.seed(16)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Please provide arguments: disease_name, PPI_file, KEGG_file, Coexpression_file, output directory.")
}

# Assign arguments to variables
disease_name <- args[1]
ppi_file <- args[2]
kegg_file <- args[3]
coexpression_file <- args[4]
output_dir <- args[5]

# Ensure the output path exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Read disease input file
disease_gene_data_frame_score <- read.table(paste0(disease_name, '_score.txt'), header = TRUE)
disease_gene_data_frame_score$P_score <- abs(disease_gene_data_frame_score$P_score)
disease_gene_data_frame_score$Lead_SNP <- gsub('-', '.', disease_gene_data_frame_score$Lead_SNP)

Layers <- list()

if (file.exists(ppi_file)) {
  t1 <- read.csv(ppi_file)
  row.names(t1) <- t1[,1]
  t1 <- t1[,2:ncol(t1)]
  t1 <- as.matrix(t1)
  PPI_cor.net <- graph.adjacency(t1, mode = "undirected", weighted = TRUE)
  PPI_cor.net <- simplify(PPI_cor.net, remove.multiple = TRUE, remove.loops = TRUE)
  V(PPI_cor.net)$type <- 'PPI'
  Layers[["PPI"]] <- PPI_cor.net
  print(paste("PPI Number of Edges:", ecount(PPI_cor.net)))
  print(paste("PPI Number of Nodes:", vcount(PPI_cor.net)))
}

if (file.exists(kegg_file)) {
  t2 <- read.csv(kegg_file)
  row.names(t2) <- t2[,1]
  t2 <- t2[,2:ncol(t2)]
  t2 <- as.matrix(t2)
  Kegg_cor.net <- graph.adjacency(t2, mode = "directed", weighted = TRUE)
  Kegg_cor.net <- simplify(Kegg_cor.net, remove.multiple = TRUE, remove.loops = TRUE)
  V(Kegg_cor.net)$type <- 'Kegg'
  Layers[["KEGG"]] <- Kegg_cor.net
  print(paste("KEGG Number of Edges:", ecount(Kegg_cor.net)))
  print(paste("KEGG Number of Nodes:", vcount(Kegg_cor.net)))
}

if (file.exists(coexpression_file)) {
  t3 <- read.csv(coexpression_file, header = TRUE)
  row.names(t3) <- t3[,1]
  t3 <- t3[,2:ncol(t3)]
  t3 <- as.matrix(t3)
  Coexpression_cor.net <- graph.adjacency(t3, mode = "undirected", weighted = TRUE)
  Coexpression_cor.net <- simplify(Coexpression_cor.net, remove.multiple = TRUE, remove.loops = TRUE)
  V(Coexpression_cor.net)$type <- 'Coexpression'
  Layers[["Coexpression"]] <- Coexpression_cor.net
  print(paste("Coexpression Number of Edges:", ecount(Coexpression_cor.net)))
  print(paste("Coexpression Number of Nodes:", vcount(Coexpression_cor.net)))
}

L <- length(Layers)
Network_List <- names(Layers)
Node_Names_all <- character()
for (i in 1:L) {
  Node_Names_Layer <- V(Layers[[i]])$name
  Node_Names_all <- c(Node_Names_all,Node_Names_Layer)
}
  
Node_Names_all <- unique(Node_Names_all)
Layers_New <- vector("list", L)
  
for (i in 1:L){
  Node_Names_Layer <- V(Layers[[i]])$name
  Missing_Nodes <- Node_Names_all [which(!Node_Names_all %in% Node_Names_Layer)]
  Layers_New[[i]] <- add_vertices(Layers[[i]] ,length(Missing_Nodes), name=Missing_Nodes)
}
  
  
vector_check <- numeric(length = L)  

for (i in 1:L){
  vector_check[i] <- vcount(Layers_New[[i]])  
}

table(all(vector_check == vector_check[1]))

N <- length(Node_Names_all)

Idem_Matrix <- Diagonal(N, x = 1)

SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L)
Col_Node_Names <- character()
Row_Node_Names <- character()
Degree <- vector("list", L)
delta <- 0.5
print('Start building network')
for (i in 1:L){
  names(Layers_New)[i] <- Network_List[i]
  Adjacency_Layer <- as_adjacency_matrix(Layers_New[[i]],attr = 'weight',sparse = TRUE)
  Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
  Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
  Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
  Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
  Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
  
  #Modification:------------------------------------------------------------------------------------------------------------------------------------------------- 
  #Calculate in-degree and out-degree of the adjacency matrix
  InDegree <- rowSums(Adjacency_Layer)
  OutDegree <- colSums(Adjacency_Layer)
  #if the layer is undirected, perform symmetric normalization on the adjacency matrix. Otherwise, perform column normalization (normalized based on end points)
  if (is.directed(Layers_New[[i]]) == FALSE) {
    D2 <- diag(1/sqrt(apply(Adjacency_Layer, 2, sum)))
    D2[!is.finite(D2)] <- 0
    Adjacency_Layer_Normalized <- D2 %*% Adjacency_Layer %*% D2
    DegreeVector <- (InDegree + OutDegree) / 2
  } else {
    D_in <- diag(1/apply(Adjacency_Layer, 2, sum))
    D_in[!is.finite(D_in)] <- 0
    Adjacency_Layer_Normalized <- Adjacency_Layer %*% D_in
    DegreeVector <- InDegree + OutDegree
  }
  Position_ini_row <- 1 + (i-1)*N
  Position_end_row <- N + (i-1)*N
  SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*(Adjacency_Layer_Normalized)
  Degree[[i]] <- DegreeVector
}
# Modification: calculate spearman correlation between layers-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
InterSpearman <- Matrix(0,ncol=L,nrow=L)
for (l1 in 1:L) {
  for (l2 in 1:L) {
    InterSpearman[l1, l2] <- cor(x = Degree[[l1]], Degree[[l2]], method = "spearman")
    InterSpearman[l2, l1] <- InterSpearman[l1, l2]
  }
}

for (i in 1:L){
  Position_ini_row <- 1 + (i-1)*N
  Position_end_row <- N + (i-1)*N
  for (j in 1:L){
    Position_ini_col <- 1 + (j-1)*N
    Position_end_col <- N + (j-1)*N
    if (j != i){
      SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix *InterSpearman[i,j]
    }
  }
}
rownames(SupraAdjacencyMatrix) <- Row_Node_Names
colnames(SupraAdjacencyMatrix) <- Col_Node_Names

tau <- rep(1,L)/L
Seeds_Genes_Scores <- numeric(length = length(disease_gene_data_frame_score$Lead_SNP)*L)
Seed_Genes_Layer_Labeled <- character(length = length(disease_gene_data_frame_score$Lead_SNP)*L)

Current_Gene_Labeled <- character()

for (k in 1:L){
  Current_Gene_Labeled <- c(Current_Gene_Labeled,paste(disease_gene_data_frame_score$Lead_SNP[k],k,sep="_",collapse = "") )
  for (j in 1:length(disease_gene_data_frame_score$Lead_SNP)){
    Seed_Genes_Layer_Labeled[((k-1)*length(disease_gene_data_frame_score$Lead_SNP))+ j] <- paste(disease_gene_data_frame_score$Lead_SNP[j],k,sep="_",collapse = "") 
    Seeds_Genes_Scores[((k-1)*length(disease_gene_data_frame_score$Lead_SNP))+ j] <- tau[k]*disease_gene_data_frame_score$P_score[j]/sum(disease_gene_data_frame_score$P_score)
  }
}  
Seeds_Score <- data.frame(Seeds_ID = Seed_Genes_Layer_Labeled,
                         Score = Seeds_Genes_Scores, stringsAsFactors = FALSE)

NetworkSize <- ncol(SupraAdjacencyMatrix)

prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)

seed_pos <- which(colnames(SupraAdjacencyMatrix) %in% Seeds_Score[,1])
for (i in seed_pos){
  prox_vector[i] <- Seeds_Score$Score[which(Seeds_Score$Seeds_ID == colnames(SupraAdjacencyMatrix)[i])]
}
prox_vector  <- prox_vector/sum(prox_vector)
restart_vector <-  prox_vector

rwr <- function(SupraAdjacencyMatrix, prox_vector, alpha, threshold, err, iter){
  restart_vector <-  prox_vector
  while(err >= threshold){
    old_prox_vector <- prox_vector
    prox_vector <- (1-alpha)*(SupraAdjacencyMatrix  %*% prox_vector) + alpha*restart_vector
    err <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1;
  }
  return (prox_vector)
}
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}  
err <- 1
iter <- 1
threshold <- 1e-15
alpha <- 0.7

prox_vector <- as.matrix(rwr(SupraAdjacencyMatrix, prox_vector, alpha, threshold, err,iter))
Random_Walk_Results <- as.data.frame(cbind(rownames(prox_vector),prox_vector))
Random_Walk_Results$GeneNames <- gsub("_[^_]*$","",Random_Walk_Results$V1)
Random_Walk_Results$log_value <- log(as.numeric(Random_Walk_Results$V2)+1)

rank_global_names <- character(length = length(unique(Random_Walk_Results$GeneNames)))
rank_global_scores <- numeric(length = length(unique(Random_Walk_Results$GeneNames)))
for ( i in (1:length(unique(Random_Walk_Results$GeneNames)))){
  
  rank_global_names[i] <- as.character(Random_Walk_Results$GeneNames[i])
  rank_global_scores[i] <- exp(mean(Random_Walk_Results$log_value[which(Random_Walk_Results$GeneNames == rank_global_names[i])]))
  
}
rank_global <- as.data.frame(cbind(rank_global_names,rank_global_scores))
rank_global$rank_global_scores <- as.numeric(rank_global$rank_global_scores)
rank_global$rank_global_scores <- minMax(rank_global$rank_global_scores)
write.table(rank_global, paste0(output_dir,'/',disease_name,'_RWR_M_result.txt'),quote = F,row.names=F,col.names=T)

## We sort the genes according to their score. 

rank_global_NoSeeds <- rank_global[which(!rank_global$rank_global_names %in% disease_gene_data_frame_score$Lead_SNP),]
write.table(rank_global_NoSeeds,paste0(output_dir,'/',disease_name,'_No_Seeds_RWR_M_result.txt'),quote = F,row.names=F,col.names=T)

print('Finish RWR. Starting permutation test')

# Permutation Test

alpha <- 0.7
nperms <- 1000
obs_rwr <- rank_global # compute observed RWR scores
perm_rwr <- matrix(0, nrow = nperms, ncol = N) # initialize permutation matrix

for (j in 1:nperms) {
  err <- 1
  iter <- 1
  threshold <- 1e-15
  
  # Randomly permute the seed nodes (uncommend to select option)
  
  # Option 1: Randomly reassign existing values for seed scores by shuffling restart_vector values.
  s_perm <- sample(restart_vector, length(restart_vector))

  # Option 2: Randomly select new seed nodes and assign uniform weights.
  # Uncomment the lines below to enable this option:
  #perm_seeds <- sample(seq_len(N), size = length(seed_pos), replace = FALSE) #
  #s_perm <- numeric(N * L)
  #s_perm[perm_seeds] <- Seeds_Score$Score / sum(Seeds_Score$Score)
  
  # Perform RWR with the permuted restart vector
  perm_rwr_ind <- as.matrix(rwr(SupraAdjacencyMatrix, s_perm, alpha, threshold, err, iter))
  perm_rwr_ind <- as.data.frame(cbind(rownames(perm_rwr_ind), perm_rwr_ind))
  
  # Process the permuted RWR results to calculate scores
  perm_rwr_ind$V1 <- gsub("_[^_]*$", "", perm_rwr_ind$V1)
  perm_rwr_ind$V2 <- log(as.numeric(perm_rwr_ind$V2) + 1)
  
  perm_rwr_ind_names <- unique(perm_rwr_ind$V1)
  perm_rwr_ind_scores <- sapply(perm_rwr_ind_names, function(gene) {
    exp(mean(perm_rwr_ind$V2[perm_rwr_ind$V1 == gene]))
  })
  
  # Normalize scores between 0 and 1
  perm_rwr_ind_scores <- minMax(perm_rwr_ind_scores)
  perm_rwr[j, ] <- perm_rwr_ind_scores # store the permuted scores
}

# Calculate p-values based on the permutation distributions
colnames(perm_rwr) <- perm_rwr_ind_names
pvals <- sapply(1:dim(obs_rwr)[1], function(i) {
  sum(perm_rwr[, i] >= obs_rwr$rank_global_scores[i]) / nperms
})


pval_result <- data.frame(rank_global_names = perm_rwr_ind_names, Pvalue = pvals)

# Merge p-values with observed RWR results and sort
rank_global_result <- merge(rank_global, pval_result, by = 'rank_global_names', all = TRUE)
rank_global_result$Pvalue <- as.numeric(rank_global_result$Pvalue)
rank_global_result_sort <- rank_global_result[order(rank_global_result$rank_global_scores, decreasing = TRUE), ]
rank_global_result_sort_no_seed <- rank_global_result_sort[!rank_global_result_sort$rank_global_names %in% disease_gene_data_frame_score$Lead_SNP, ]

write.table(rank_global_result_sort, paste0(output_dir, '/', disease_name, '_permutation_result.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(rank_global_result_sort_no_seed, paste0(output_dir, '/', disease_name, '_permutation_no_seed_result.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE)

print(paste0(disease_name, '_RWR finished!'))
