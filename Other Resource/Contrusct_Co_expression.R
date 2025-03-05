library(reshape)
library(igraph)
library("biomaRt")
data.tissue <-read.csv('RNA_consensus_tissue_gene_data/rna_tissue_consensus.tsv',head=T,sep='\t')
data.celline = read.csv('RNA_HPA_immune_cell_sample_gene_data/rna_immune_cell_sample.tsv',head=T,sep='\t')
data.celline = data.celline[c(4,5,3,8)]
colnames(data.tissue) = c('Ensemble_ID', 'Gene.name', 'Sample', 'Value')
colnames(data.celline) = c('Ensemble_ID', 'Gene.name', 'Sample', 'Value')
data.tissue = data.tissue[which((data.tissue$Sample=='lung')|(data.tissue$Sample=='kidney')),]
data.celline=data.celline[which((data.celline$Sample=='NK-cell')|(data.celline$Sample=='T-reg')|(data.celline$Sample=='naive B-cell')|(data.celline$Sample=='naive CD8 T-cell')|(data.celline$Sample=='memory CD8 T-cell')|(data.celline$Sample=='naive CD4 T-cell')|(data.celline$Sample=='memory CD4 T-cell')),]
data_Coexpr <- rbind(data.tissue,data.celline)

check.HGNC.symbols <- function(GeneNames){
  
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  ensemble_HGNC <- getBM(attributes=c('hgnc_symbol'), 
                         filters='hgnc_symbol', values=GeneNames, mart=ensembl)
  
  return(ensemble_HGNC$hgnc_symbol)
}
ensemble_HGNC <- check.HGNC.symbols(unique(data_Coexpr$Gene.name))
data3_Coexpr <- data_Coexpr[which(data_Coexpr$Gene.name %in% ensemble_HGNC),]

matrix_Coexpr<-cast(data3_Coexpr, Gene.name ~ Sample,mean,value="Value") 
rownames(matrix_Coexpr)<-matrix_Coexpr[,1]
matrix_Coexpr[,1]<-NULL
#remove rows with NA mean
c = which(is.na(matrix_Coexpr))
matrix_Coexpr = matrix_Coexpr[which(!rownames(matrix_Coexpr)%in% c),]

# Identify rows (or columns) with zero standard deviation
sd_values <- apply(matrix_Coexpr, 1, sd)  # For rows
zero_sd_rows <- which(sd_values == 0)

sd_values_col <- apply(matrix_Coexpr, 2, sd)  # For columns
zero_sd_cols <- which(sd_values_col == 0)

# Display problematic rows and columns
print("Rows with zero standard deviation:")
print(zero_sd_rows)

print("Columns with zero standard deviation:")
print(zero_sd_cols)

# Remove rows with zero standard deviation
if (length(zero_sd_rows) > 0) {
  matrix_Coexpr <- matrix_Coexpr[-zero_sd_rows, ]
}

# Remove columns with zero standard deviation
if (length(zero_sd_cols) > 0) {
  matrix_Coexpr <- matrix_Coexpr[, -zero_sd_cols]
}

# Transpose if needed and compute the correlation matrix
CR <- cor(t(matrix_Coexpr), method = "spearman")

diag(CR) = 0
inds <- which(abs(CR[1:nrow(CR),]) >= 0.75, arr.ind=TRUE)
inds = as.data.frame(inds)
cor_value = c()
for (i in 1:dim(inds)[1]){
  cor_value[i] = CR[inds$row[i],inds$col[i]]
}
inds$weight = abs(cor_value)
row_name =c()
col_name = c()

for (i in 1:dim(inds)[1]){
  row_name[i] = rownames(CR)[inds$row[i]]
  col_name[i] = colnames(CR)[inds$col[i]]
}
inds$row_name = row_name
inds$col_name = col_name

inds$network_weight=(inds$weight/min(inds$weight))
co_express_network_edge_list = inds[c(5,6,4)]
write.csv(co_express_network_edge_list,'Coexpression_AAV_edgelist_corr0.75_with_weight.csv',quote=F,row.names=F,col.names=T)

