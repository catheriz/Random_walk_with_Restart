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
set.seed(16)
disease_name_list = c('PM','DM','JDM','Jo1','total_IIMs')
for (disease_name in disease_name_list){
  #Read disease input file
  disease_gene_data_frame_score = read.table(paste0(disease_name,'_score.txt'),head=T)
  disease_gene_data_frame_score$P_score = abs(disease_gene_data_frame_score$P_score)
  disease_gene_data_frame_score$Lead_SNP = gsub('-','.',disease_gene_data_frame_score$Lead_SNP )
  
  #Load PPI network
  t1=read.csv(paste0(disease_name,'_PPI_adjacency_matrix.csv'))
  row.names(t1) = t1[,1]
  t1=t1[,2:dim(t1)[2]]
  t1 = as.matrix(t1)
  #isSymmetric(t1, check.attributes = FALSE)
  PPI_cor.net = graph.adjacency(t1, mode="undirected", weighted = TRUE)
  PPI_cor.net= simplify(PPI_cor.net,remove.multiple = TRUE, remove.loops = TRUE)
  V(PPI_cor.net)$type = 'PPI' 
  print(paste("PPI Number of Edges: ", ecount(PPI_cor.net),sep="",collapse = NULL))
  print(paste("PPI Number of Nodes: ", vcount(PPI_cor.net),sep="",collapse = NULL))
  #Load Kegg Network
  t2 = read.csv(paste0(disease_name,'_Kegg_gene_adjacency_matrix_directed.csv'))
  row.names(t2) = t2[,1]
  t2=t2[,2:dim(t2)[2]]
  t2 = as.matrix(t2)
  #isSymmetric(t2, check.attributes = FALSE)
  Kegg_cor.net = graph.adjacency(t2, mode="directed", weighted = TRUE)
  Kegg_cor.net= simplify(Kegg_cor.net,remove.multiple = TRUE, remove.loops = TRUE)
  print(paste("KEGG Number of Edges: ", ecount(Kegg_cor.net),sep="",collapse = NULL))
  print(paste("KEGG Number of Nodes: ", vcount(Kegg_cor.net),sep="",collapse = NULL))
  V(Kegg_cor.net)$type = 'Kegg'
  
  #Load Co-expression Network
  t3 = read.csv('Coexpression_immune_tissue_cell_line_network_edgelist_corr0.75_adjacency_matrix_undirected.csv',head=T)
  rownames(t3)=t3[,1]
  t3 = t3[c(2:(dim(t3)[2]))]
  t3 = as.matrix(t3)
  Coexpression_cor.net = graph.adjacency(t3, mode="undirected",weighted = TRUE)
  Coexpression_cor.net = simplify(Coexpression_cor.net,remove.multiple = TRUE, remove.loops = TRUE)
  print(paste("Coexpression Number of Edges: ", ecount(Coexpression_cor.net),sep="",collapse = NULL))
  print(paste("Coexpression Number of Nodes: ", vcount(Coexpression_cor.net),sep="", collapse = NULL))
  V(Coexpression_cor.net)$type = 'Coexpression'
  
  # The following code was adopted from RWR-M source code
  Network_List = c("PPI","KEGG","Coexpression")
  L = length(Network_List)
  
  Layers = vector("list", L)
  
  for (i in 1:L){
    if (Network_List[i]=="PPI"){
      Layers[[i]] = PPI_cor.net
      names(Layers)[i] = "PPI"
    } else {

      if (Network_List[i]=="KEGG"){  
        Layers[[i]] = Kegg_cor.net
        names(Layers)[i] = "KEGG"
        
      } else {
        
        if (Network_List[i]=="Coexpression"){
          Layers[[i]] = Coexpression_cor.net
          names(Layers)[i] = "Coexpression"
     
        } 
      } 
    }
  }
  
  Node_Names_all = character()
  for (i in 1:L) {
    Node_Names_Layer = V(Layers[[i]])$name
    Node_Names_all =c(Node_Names_all,Node_Names_Layer)
  }
  
  Node_Names_all = unique(Node_Names_all)
  Layers_New = vector("list", L)
  
  for (i in 1:L){
    Node_Names_Layer = V(Layers[[i]])$name
    Missing_Nodes = Node_Names_all [which(!Node_Names_all %in% Node_Names_Layer)]
    Layers_New[[i]] = add_vertices(Layers[[i]] ,length(Missing_Nodes), name=Missing_Nodes)
  }
  
  vector_check = numeric(length = L)  
  
  for (i in 1:L){
    vector_check[i] = vcount(Layers_New[[i]])  
  }
  
  table(all(vector_check == vector_check[1]))
  
  N = length(Node_Names_all)

  Idem_Matrix = Diagonal(N, x = 1)
  
  SupraAdjacencyMatrix = Matrix(0,ncol=N*L,nrow=N*L)
  Col_Node_Names = character()
  Row_Node_Names = character()
  Degree = vector("list", L)
  delta = 0.5
  
  for (i in 1:L){
    names(Layers_New)[i] = Network_List[i]
    Adjacency_Layer =  as_adjacency_matrix(Layers_New[[i]],attr = 'weight',sparse = TRUE)
    Adjacency_Layer = Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Col_Names = paste(colnames(Adjacency_Layer),i,sep="_")
    Layer_Row_Names = paste(rownames(Adjacency_Layer),i,sep="_")
    Col_Node_Names = c(Col_Node_Names,Layer_Col_Names)
    Row_Node_Names = c(Row_Node_Names,Layer_Row_Names)
    
    #Modification:------------------------------------------------------------------------------------------------------------------------------------------------- 
    #Calculate in-degree and out-degree of the adjacency matrix
    InDegree = rowSums(Adjacency_Layer)
    OutDegree = colSums(Adjacency_Layer)
    #if the layer is undirected, perform symmetric normalization on the adjacency matrix. Otherwise, perform column normalization (normalized based on end points)
    if (is.directed(Layers_New[[i]]) == FALSE) {
      D2 = diag(1/sqrt(apply(Adjacency_Layer, 2, sum)))
      D2[!is.finite(D2)] = 0
      Adjacency_Layer_Normalized= D2 %*% Adjacency_Layer %*% D2
      DegreeVector = (InDegree + OutDegree) / 2
    } else {
      D_in= diag(1/apply(Adjacency_Layer, 2, sum))
      D_in[!is.finite(D_in)] = 0
      Adjacency_Layer_Normalized = Adjacency_Layer %*% D_in
      DegreeVector = InDegree + OutDegree
    }
    Position_ini_row = 1 + (i-1)*N
    Position_end_row = N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] = (1-delta)*(Adjacency_Layer_Normalized)
    Degree[[i]] = DegreeVector
  }
  # Modification: calculate spearman correlation between layers-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  InterSpearman= Matrix(0,ncol=L,nrow=L)
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      InterSpearman[l1, l2] =cor(x = Degree[[l1]], Degree[[l2]], method = "spearman")
      InterSpearman[l2, l1] = InterSpearman[l1, l2]
    }
  }
  
  for (i in 1:L){
    Position_ini_row = 1 + (i-1)*N
    Position_end_row = N + (i-1)*N
    for (j in 1:L){
      Position_ini_col = 1 + (j-1)*N
      Position_end_col = N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] = (delta/(L-1))*Idem_Matrix *InterSpearman[i,j]
      }
    }
  }
  rownames(SupraAdjacencyMatrix) = Row_Node_Names
  colnames(SupraAdjacencyMatrix) = Col_Node_Names
  
  tau = rep(1,L)/L
  Seeds_Genes_Scores = numeric(length = length(disease_gene_data_frame_score$Lead_SNP)*L)
  Seed_Genes_Layer_Labeled = character(length = length(disease_gene_data_frame_score$Lead_SNP)*L)
  
  Current_Gene_Labeled = character()
  
  for (k in 1:L){
    Current_Gene_Labeled = c(Current_Gene_Labeled,paste(disease_gene_data_frame_score$Lead_SNP[k],k,sep="_",collapse = "") )
    for (j in 1:length(disease_gene_data_frame_score$Lead_SNP)){
      Seed_Genes_Layer_Labeled[((k-1)*length(disease_gene_data_frame_score$Lead_SNP))+ j] = paste(disease_gene_data_frame_score$Lead_SNP[j],k,sep="_",collapse = "") 
      Seeds_Genes_Scores[((k-1)*length(disease_gene_data_frame_score$Lead_SNP))+ j] = tau[k]*disease_gene_data_frame_score$P_score[j]/sum(disease_gene_data_frame_score$P_score)
    }
  }  
  Seeds_Score = data.frame(Seeds_ID = Seed_Genes_Layer_Labeled,
                            Score = Seeds_Genes_Scores, stringsAsFactors = FALSE)
  
  NetworkSize = ncol(SupraAdjacencyMatrix)
  
  prox_vector = matrix(0,nrow = NetworkSize,ncol=1)
  
  seed_pos = which(colnames(SupraAdjacencyMatrix) %in% Seeds_Score[,1])
  for (i in seed_pos){
    prox_vector[i] = Seeds_Score$Score[which(Seeds_Score$Seeds_ID == colnames(SupraAdjacencyMatrix)[i])]
  }
  prox_vector  = prox_vector/sum(prox_vector)
  restart_vector =  prox_vector
  
  rwr <- function(Adjacency_matrix, prox_vector, alpha, Threeshold, err,iter){
    restart_vector =  prox_vector
    while(err >= Threeshold){
      old_prox_vector = prox_vector
      prox_vector = (1-alpha)*(SupraAdjacencyMatrix  %*% prox_vector) + alpha*restart_vector
      err = sqrt(sum((prox_vector-old_prox_vector)^2))
      iter = iter + 1;
    }
    return (prox_vector)
  }
  minMax <- function(x) {
   (x - min(x)) / (max(x) - min(x))
  }  
  err = 1
  iter = 1
  Threeshold = 1e-15
  alpha = 0.7
  
  prox_vector = as.matrix(rwr(SupraAdjacencyMatrix, prox_vector, alpha, Threeshold, err,iter))
  Random_Walk_Results = as.data.frame(cbind(rownames(prox_vector),prox_vector))
  Random_Walk_Results$GeneNames = gsub("_[^_]*$","",Random_Walk_Results$V1)
  Random_Walk_Results$log_value = log(as.numeric(Random_Walk_Results$V2)+1)
  
  rank_global_names = character(length = length(unique(Random_Walk_Results$GeneNames)))
  rank_global_scores = numeric(length = length(unique(Random_Walk_Results$GeneNames)))
  for ( i in (1:length(unique(Random_Walk_Results$GeneNames)))){
    
    rank_global_names[i] = as.character(Random_Walk_Results$GeneNames[i])
    rank_global_scores[i] = exp(mean(Random_Walk_Results$log_value[which(Random_Walk_Results$GeneNames == rank_global_names[i])]))
    
  }
  rank_global = as.data.frame(cbind(rank_global_names,rank_global_scores))
  rank_global$rank_global_scores = as.numeric(rank_global$rank_global_scores)
  rank_global$rank_global_scores = minMax(rank_global$rank_global_scores)
  write.table(rank_global, paste0('RWR_Zscore_norm_modification_no_final_col/',disease_name,'_RWR_M_method_3_layer_PPI_KEGG_Coexpression_gene_1e-15.txt'),quote = F,row.names=F,col.names=T)
  
  ## We sort the genes according to their score. 
  
  rank_global_NoSeeds = rank_global[which(!rank_global$rank_global_names %in% disease_gene_data_frame_score$Lead_SNP),]
  write.table(rank_global_NoSeeds,paste0('RWR_Zscore_norm_modification_no_final_col/',disease_name,'_No_Seeds_RWR_M_method_3_layer_PPI_KEGG_Coexpression_gene_1e-15.txt'),quote = F,row.names=F,col.names=T)
  
  
  
  #Permutation Test
  
  alpha = 0.7
  nperms = 1000
  obs_rwr = rank_global # compute observed RWR scores
  perm_rwr =matrix(0, nrow=nperms, ncol=N) # initialize permutation matrix
  for (j in 1:nperms) {
    err = 1
    iter = 1
    Threeshold = 1e-15
    s_perm <- sample(restart_vector, length(restart_vector)) # randomly permute the seed nodes
    perm_rwr_ind = as.matrix(rwr(SupraAdjacencyMatrix, s_perm, alpha, Threeshold, err,iter))
    perm_rwr_ind = as.data.frame(cbind(rownames(perm_rwr_ind),perm_rwr_ind))
    perm_rwr_ind $V1 = gsub("_[^_]*$","",Random_Walk_Results$V1)
    perm_rwr_ind $V2 = log(as.numeric(perm_rwr_ind$V2)+1)
    perm_rwr_ind_names = character(length = length(unique(perm_rwr_ind $V1)))
    perm_rwr_ind_scores = numeric(length = length(unique(perm_rwr_ind $V1)))
    for ( i in (1:length(unique(perm_rwr_ind $V1)))){
      perm_rwr_ind_names[i] = as.character(perm_rwr_ind $V1[i])
      perm_rwr_ind_scores[i] = exp(mean(perm_rwr_ind $V2[which(perm_rwr_ind $V1 == perm_rwr_ind_names[i])]))
    }
    perm_rwr_ind_scores = minMax(perm_rwr_ind_scores)
    perm_rwr[j,] =perm_rwr_ind_scores # compute RWR scores on permuted seed set
  }
  #Calculate the geometric term 
  colnames(perm_rwr) = perm_rwr_ind_names
  pvals = list()
  for (i in 1:dim(obs_rwr)[1]){
    num_sum = sum(perm_rwr[,i]>= obs_rwr$rank_global_scores[i])
    pvals[i] = as.numeric(num_sum/nperms)
    }
  
  pval_result = as.data.frame(cbind(perm_rwr_ind_names,pvals))
  colnames(pval_result) = c('rank_global_names','Pvalue')
  
  rank_global_result = merge(rank_global,pval_result,by='rank_global_names',all=T)
  rank_global_result$Pvalue = as.numeric(unlist(rank_global_result$Pvalue))
  rank_global_result_sort = rank_global_result[order(rank_global_result$rank_global_scores,decreasing = TRUE),]
  rank_global_result_sort_no_seed = rank_global_result_sort[which(!rank_global_result_sort$rank_global_names %in% disease_gene_data_frame_score$Lead_SNP),]
  
  write.table(rank_global_result_sort,paste0('RWR_Zscore_norm_modification_no_final_col/',disease_name,'_permutation_result.txt'),quote = F,row.names=F,col.names=F)
  
  write.table(rank_global_result_sort_no_seed,paste0('RWR_Zscore_norm_modification_no_final_col/',disease_name,'_permutation_no_seed_result.txt'),quote = F,row.names=F,col.names=F)
  rank_global_result_sort_no_seed_sig = rank_global_result_sort_no_seed[which(rank_global_result_sort_no_seed$Pvalue<=0.05),]
  rank_global_result_sort_no_seed_sig$Pvalue = as.numeric(rank_global_result_sort_no_seed_sig$Pvalue)
  
  write.table(rank_global_result_sort_no_seed_sig,paste0('RWR_Zscore_norm_modification_no_final_col/',disease_name,'_permutation_result_significant.txt'),quote = F,row.names=F,col.names=F)
  c = rank_global_result_sort_no_seed_sig$rank_global_names[1:30]
  c = append(c,disease_gene_data_frame_score$Lead_SNP)
  write.table(c, paste0('RWR_Zscore_norm_modification_no_final_col/Top30_genes_and_disease_genes_',disease_name,'_3_Network_PPI_KEGG_Coexpression_immune_gene_1e-15.txt'),quote = F,row.names=F,col.names=F)
}

