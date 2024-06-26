# combine and normalize Osteo/Adipo and Chondro lineage

# data sets pre-filtered (remove high percent MT, low feature/read count)

# RNA data normalization procedure:
# 1. Begin with raw RNA read counts
# 2. Merge counts from each dataset (O/A, C), only including shared genes measured in all datasets
# 3. Split object by Batch.ID
# 4. Perform SCTransform normalization on each batch separately
# 5. Re-merge SCT data, again only with common genes reported by SCT
# 6. Perform quantile normalization on SCT data
# 7. Lastly, score cell cycle and keep only G1 cells


library(tidyverse)
library(Seurat)


# load Osteo/Adipo MTL data
MTL_obj = read_rds("Data/MTL/2020-06-25 MTLv4 integrated harmony.RDS")
dim(MTL_obj) # 22194 genes x 31527 cells (SCT slot)


# load chondrocyte lineage data (GSE160625)
chondro_obj = readRDS("Data/GSE160625_Chondro/GSE160625_filtnorm_C59subset_seurat.RDS")
chondro_obj$orig.ident = str_replace(chondro_obj$orig.ident, "GSM.+_(?=D)","Chondro ") # rename GSE160625 labels to Chondro
chondro_obj$orig.ident = str_replace(chondro_obj$orig.ident, "GSM.+_Cp","Chondro D0") # Cp is day 0 of chondro differentiation
dim(chondro_obj) # 24506 genes x 11136 cells (SCT slot)



# add source metadata before merging
MTL_obj$source = "MTL"
chondro_obj$source = "chondro"

# set all active assay to RNA (so we merge the raw RNA counts) and remove SCT assay for memory (will recompute later)
if (MTL_obj@active.assay == "SCT") {MTL_obj@active.assay = "RNA"; MTL_obj[['SCT']] <- NULL}
if (chondro_obj@active.assay == "SCT") {chondro_obj@active.assay = "RNA"; chondro_obj[['SCT']] <- NULL}
dim(MTL_obj) # 33538 genes x 31527 cells (RNA slot)
dim(chondro_obj) # 45924 genes x 11136 cells (RNA slot)

# merge into one seurat object
obj_merge = merge(MTL_obj, chondro_obj)
rm(MTL_obj, chondro_obj); gc() # to free some memory
dim(obj_merge) # 57689 genes x 42663 cells


# update some metadata (Lineage and Batch.ID)
obj_merge$Lineage[obj_merge$source == "chondro"] = "Chondro"

# assign the chondro and PDX data distinct Batch.IDs
obj_merge$Batch.ID[obj_merge$source == "chondro"] = 4 # assign chondro data to a new batch id 4

# split data to run SCTransform on each batch
obj.list <- SplitObject(obj_merge, split.by = "Batch.ID")
rm(obj_merge); gc() # to free some memory
obj.list <- lapply(obj.list, function(x) { x <- SCTransform(x, return.only.var.genes = F, vst.flavor="v2") }) # using SCTransform v2


# merge SCT data back together
obj_merge_sct = merge(obj.list[[1]], obj.list[-1])
VariableFeatures(obj_merge_sct) <- Reduce(intersect, lapply(obj.list, VariableFeatures))
rm(obj.list); gc() # to free some memory
dim(obj_merge_sct) # 28007 genes x 42663 cells


# remove some columns from metadata (not relevant to all datasets)
rem_cols = c("gemgroup", "old.ident", "discard.mito", "discard.features", "discard.counts",
             "CC.Difference", "SCT_snn_res.0.3", "SCT_snn_res.0.1", "SCT_snn_res.0.05", "SCT_snn_res.0.02", "seurat_clusters")
for (ic in rem_cols) { obj_merge_sct[[ic]] <- NULL }

# cell cycle scoring and filtering
obj_merge_sct = CellCycleScoring(obj_merge_sct,
                                 s.features = union(cc.genes$s.genes, cc.genes.updated.2019$s.genes),
                                 g2m.features = union(cc.genes$g2m.genes, cc.genes.updated.2019$g2m.genes) )


# save data
savedir = "Data/MTL_OAC_SCT_qnorm"
saveRDS(obj_merge_sct, file.path(savedir,"MTL_OAC_merged_SCT.rds"))
#obj_merge_sct = readRDS( file.path(savedir,"MTL_OAC_merged_SCT.rds") ) # to load



# save data for NMF analysis
savedir = "Data/MTL_OAC_SCT_qnorm"
obj_merge_sct = readRDS( file.path(savedir,"MTL_OAC_merged_SCT.rds") ) # to load
meta_all = obj_merge_sct@meta.data

# quantile normalize
library(preprocessCore)
d = obj_merge_sct@assays$SCT@data # pull SCT data slot
rm(obj_merge_sct); gc() # remove Seurat object to save memory
dq = normalize.quantiles(d %>% as.matrix, keep.names = T)

# subset G1 phase cells
keep_phase = meta_all$Phase == "G1"
dq_filt = dq[,keep_phase]
meta_filt = meta_all[keep_phase,]


# function to save all data as .csv into multiple files separated by condition
save_each_condition <- function(expr, meta, expr_dir, prec=5) {  # prec is precision for tsv (number of digits, default 5)
  if (!dir.exists(expr_dir)) { dir.create(expr_dir) } # create directory if it does not exist
  # loop through each condition in meta and save separately
  for (cond in unique(meta$orig.ident)) {
    inds = which(meta$orig.ident == cond); print(paste(cond,":",length(inds),"cells"))
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
save_each_condition(dq, meta_all, file.path(savedir,"MTL_OAC_expr_condition"))
save_each_condition(dq_filt, meta_filt, file.path(savedir,"MTL_OAC_expr_condition_ccfilt"))

