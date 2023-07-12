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
#obj_merge_sct = readRDS( file.path(savedir,"OSPDX_SCT.rds") ) # to load

# filter cell cycle phase to select only cells in G1 phase
obj_merge_filt = subset(obj_merge_sct, subset = Phase == "G1")
dim(obj_merge_filt) # 24756 genes x 12016 cells
saveRDS(obj_merge_filt, file.path(savedir,"OSPDX_SCT_ccfilt.rds"))
#obj_merge_filt = readRDS( file.path(savedir,"OSPDX_merged_SCT_ccfilt.rds") ) # to load


# quantile normalize (using target distribution)
library(preprocessCore)

# define target distribution from MTL dataset
savedir_target = "Data/MTL_OAC_SCT_qnorm"
MTL_obj = readRDS( file.path(savedir_target,"MTL_OAC_merged_SCT.rds") )
common.genes = intersect(rownames(obj_merge_sct), rownames(MTL_obj)) # subset only the genes from both MTL and PDX datasets
length(common.genes) # 20201 genes
d_target = MTL_obj@assays$SCT@data[common.genes,] # pull SCT data slot for expression target distributions
target = normalize.quantiles.determine.target(d_target %>% as.matrix)

# qnorm with target distribution
d = obj_merge_sct@assays$SCT@data[common.genes,] # pull SCT data slot
dq = normalize.quantiles.use.target(d %>% as.matrix, target)
rownames(dq) = rownames(d)
colnames(dq) = colnames(d)
dim(dq) # 20201 genes x 19538 cells


# save as .csv into multiple files separated by sample
prec = 5 # precision for tsv (number of digits)
for (pdx_i in unique(obj_merge_sct$sample_id)) {
  inds = which(obj_merge_sct$sample_id == pdx_i)
  temp = dq[,inds] %>% as.matrix
  temp = floor( temp*(10^prec) )/ (10^prec) # this truncates to prec digits for csv
  filename = paste0("expm-full-",pdx_i,".txt")
  write_tsv(as.data.frame(temp),file.path(savedir,filename), col_names=T)
}
write_tsv(data.frame(SYMBOL=rownames(dq)), file.path(savedir,"gene-conversion.txt"))

# lastly save the metadata
write_tsv(obj_merge_sct@meta.data %>% rownames_to_column("barcode"), file.path(savedir,"cell_metadata.txt"))

