# prepare gene expression data for NMF
# TARGET OS dataset

# processing:
# for bulk RNA: provided as log2(CPM+1)
# quantile normalize


library(tidyverse)

# TARGET Osteosarcoma data, accessed through Xena Browser
targetos = read_delim("Data/TARGET_OS/from_Xena/TARGET-OS.star_counts.tsv")
targetos = targetos[-c(1:4),] # first 4 rows are mapping statistics

# match ensembl id to HGNC
library(biomaRt)
ensembl_ids = sub("[.][0-9]*","",targetos$Ensembl_ID)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = ensembl_ids, mart= mart)
match_genes = left_join(data.frame(ens_id = ensembl_ids),gene_IDs,
                        by = c("ens_id"="ensembl_gene_id"), multiple="first")

exp_mat = as.matrix(targetos[,-1])
rownames(exp_mat) = match_genes$hgnc_symbol


# quantile normalize (using target distribution)
library(preprocessCore)
library(Seurat)

# first load MTL dataset to define target distribution
savedir_target = "Data/MTL_OAC_SCT_qnorm"
MTL_obj = readRDS( file.path(savedir_target,"MTL_OAC_merged_SCT.rds") )
gene_names = intersect(rownames(exp_mat), rownames(MTL_obj)) # subset only the genes from both MTL and PDX datasets
length(gene_names) # 18159 genes

# define target distribution based on MTL expression of shared genes
d_target = MTL_obj@assays$SCT@data[gene_names,] # pull SCT data slot for expression target distributions
rm(MTL_obj); gc() # remove Seurat object to save memory
target = normalize.quantiles.determine.target(d_target %>% as.matrix)

# qnorm with target distribution
d = exp_mat[gene_names,]
dq = normalize.quantiles.use.target(d %>% as.matrix, target)
rownames(dq) = rownames(d)
colnames(dq) = colnames(d)
dim(dq) # 18159 genes x 88 patients

# save as .csv into multiple files separated by sample
expr_dir = "Data/TARGET_OS/TARGETOS_expr_qnorm"
if (!dir.exists(expr_dir)) { dir.create(expr_dir) } # create directory if it does not exist
prec = 5 # precision for tsv (number of digits)
temp = dq %>% as.matrix
temp = floor( temp*(10^prec) )/ (10^prec) # this truncates to prec digits for csv
write_tsv(as.data.frame(temp),file.path(expr_dir,"target_os.txt"), col_names=T)
write_tsv(data.frame(SYMBOL=rownames(dq)), file.path(expr_dir,"gene-conversion.txt"))

