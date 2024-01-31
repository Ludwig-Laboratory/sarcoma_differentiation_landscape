# SCT + quantile normalize
# 3 OS PDX tumors

# data sets pre-filtered (remove outlier cells with high percent MT, low feature/read count)
# when loading, we further subset only cells from the "Warm" dissociation protocol

# RNA data normalization procedure:
# 1. Begin with raw RNA read counts
# 2. Split object by tumor_id
# 3. Perform SCTransform normalization on each batch separately
# 4. Re-merge SCT data, taking union of genes reported by SCT
# 5. Perform quantile normalization on SCT data (with target distribution from MTL dataset)
# 6. Lastly, score cell cycle and keep only G1 cells

library(tidyverse)
library(Seurat)


# load OS PDX data: 3 osteosarcoma patient-derived xenografts (OS PDX)
OSPDX_obj <- readRDS('Data/Osteosarcoma_PDX/RDS/2021-01-05 OS_PDX_data.subset.rds')
OSPDX_obj <- subset(OSPDX_obj, sample_type=="Warm") # filter out only Warm dissociated single-cells
dim(OSPDX_obj) # 36601 genes x 19538 cells


# split data to run SCTransform on each batch
obj.list <- SplitObject(OSPDX_obj, split.by = "sample_id")
rm(OSPDX_obj); gc() # to free some memory
obj.list <- lapply(obj.list, function(x) { x <- SCTransform(x, return.only.var.genes = F, vst.flavor="v2") }) # using SCTransform v2

# merge SCT data back together
obj_merge_sct = merge(obj.list[[1]], obj.list[-1])
VariableFeatures(obj_merge_sct) <- Reduce(intersect, lapply(obj.list, VariableFeatures))
rm(obj.list); gc() # to free some memory
dim(obj_merge_sct) # 24756 genes x 19538 cells


# cell cycle scoring and filtering
obj_merge_sct = CellCycleScoring(obj_merge_sct,
                                 s.features = union(cc.genes$s.genes, cc.genes.updated.2019$s.genes),
                                 g2m.features = union(cc.genes$g2m.genes, cc.genes.updated.2019$g2m.genes) )


# save data
savedir = "Data/OSPDX_SCT_qnorm"
saveRDS(obj_merge_sct, file.path(savedir,"OSPDX_SCT.rds"))
obj_merge_sct = readRDS( file.path(savedir,"OSPDX_SCT.rds") ) # to load
dim(obj_merge_sct) # 24756 genes x 19538 cells


# quantile normalize (using target distribution)
library(preprocessCore)

# first load MTL dataset to define target distribution
savedir_target = "Data/MTL_OAC_SCT_qnorm"
MTL_obj = readRDS( file.path(savedir_target,"MTL_OAC_merged_SCT.rds") )
gene_names = intersect(rownames(obj_merge_sct), rownames(MTL_obj)) # subset only the genes from both MTL and PDX datasets
length(gene_names) # 20201 genes

# define target distribution based on MTL expression of shared genes
d_target = MTL_obj@assays$SCT@data[gene_names,] # pull SCT data slot for expression target distributions
rm(MTL_obj); gc() # remove Seurat object to save memory
target = normalize.quantiles.determine.target(d_target %>% as.matrix)

# qnorm with target distribution
d = obj_merge_sct@assays$SCT@data[gene_names,] # pull SCT data slot
dq = normalize.quantiles.use.target(d %>% as.matrix, target)
rownames(dq) = rownames(d)
colnames(dq) = colnames(d)
dim(dq) # 20201 genes x 19538 cells


# function to save all data as .csv into multiple files separated by condition
save_each_condition <- function(expr, meta, expr_dir, prec=5) {  # prec is precision for tsv (number of digits, default 5)
  if (!dir.exists(expr_dir)) { dir.create(expr_dir) } # create directory if it does not exist
  # loop through each condition in meta and save separately
  for (cond in unique(meta$lab_id)) {
    inds = which(meta$lab_id == cond); print(paste(cond,":",length(inds),"cells"))
    temp = expr[,inds] %>% as.matrix
    temp = floor( temp*(10^prec) )/ (10^prec) # this truncates to prec digits for saving in .csv format
    filename = paste0("expm-full-",cond %>% str_replace(" ",""),".txt")
    write_tsv(as.data.frame(temp), file.path(expr_dir,filename), col_names=T)
  }
  # save gene names
  write_tsv(data.frame(SYMBOL=rownames(expr)), file.path(expr_dir,"gene-conversion.txt"))
  
  # lastly save the metadata
  write_tsv(meta %>% rownames_to_column("barcode"), file.path(expr_dir,"cell_metadata.txt"))
}
save_each_condition(dq, obj_merge_sct@meta.data, file.path(savedir,"OSPDX_expr_qnorm"))

