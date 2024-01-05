# Prep data from GSE152048 into a single Seurat object
# GSE152048 citation: Zhou et al. Nat Commun 2020
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152048

library(tidyverse)
library(Seurat)

#setwd("GSE152048") # specify directory with raw data downloaded from GEO


##### LOAD AND COMBINE #####

# get unique experiment folder names (of separate 10X datasets)
files = list.files()
data_names = files %>% str_extract("BC.+") %>% na.omit %>% unique

# load each sample with Seurat Read10X
obj_list = list()
for (i in 1:length(data_names)) {
  dir_i = data_names[i]
  data_i = Read10X(dir_i)
  obj_i = CreateSeuratObject(data_i, project=data_names[i])
  obj_list[[i]] = obj_i
}

# merge into one Seurat object
obj = merge(obj_list[[1]], obj_list[-1], add.cell.ids = data_names, project="GSE152048")

# save as .RDS
saveRDS(obj, "GSE152048_merged_seurat.RDS")


##### PROCESSING #####

# pre-processing and filtering
obj = readRDS("GSE152048_merged_seurat.RDS")

# steps:
# remove QC cells
# SCTransform, UMAP
# preliminary clustering/labeling (like their paper)
# save

# set metadata
obj$tumor_id = obj$orig.ident

# perform QC filtering of low quality cells
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt", assay = 'RNA')

library(SingleCellExperiment)
library(scater)
obj.sce = as.SingleCellExperiment(obj)
discard.mito <- isOutlier(obj.sce$percent.mt, type="higher", batch=obj.sce$orig.ident)
discard.counts <- isOutlier(obj.sce$nCount_RNA, log = TRUE, type="lower", batch=obj.sce$orig.ident)
discard.features <- isOutlier(obj.sce$nFeature_RNA, log = TRUE, type="lower", batch=obj.sce$orig.ident)

# filter based on: remove outliers with >3 mad in mito-mapped reads, <3 mad in nCount (total reads) or nFeature (unique genes mapped)
# also include a hard cutoff filter of mito<25%, nCount>1000, nFeature>500
obj.filt <- subset(obj, subset = !(discard.mito | discard.counts | discard.features) &
                       percent.mt < 25 & nCount_RNA > 1000 & nFeature_RNA > 500)

# some QC plots
plotColData(obj.sce, x="orig.ident", y="percent.mt", colour_by=I(discard.mito)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + labs(fill = 'Outlier') + ylab('% Mitochondrial Mapping')
plotColData(obj.sce, x="orig.ident", y="nCount_RNA", colour_by=I(discard.counts)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + ylab('Counts') + scale_y_log10()
plotColData(obj.sce, x="orig.ident", y="nFeature_RNA", colour_by=I(discard.features)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + ylab('Genes') + scale_y_log10()
RidgePlot(obj, group.by = "orig.ident", features = c("nCount_RNA")) + ggtitle("Counts")
RidgePlot(obj, group.by = "orig.ident", features = c("percent.mt")) + ggtitle("% Mitochondrial")
FeatureScatter(obj, 'nCount_RNA', 'nFeature_RNA', pt.size = 0.1) + NoLegend() +
  xlab('% Mitochondrial Mapping') + ylab("Counts per cell") + theme(axis.text.y=element_text(size=6)) +
  #xlim(0,100) +
  scale_x_log10(limits=c(10,200000)) +
  scale_y_log10(limits=c(10,10500))


# save data
obj <- obj.filt
rm(obj.filt, obj.sce); gc() # clear variables to conserve RAM

# cell cycle scoring
s.genes <- union(cc.genes$s.genes, cc.genes.updated.2019$s.genes) # cc.genes or cc.genes.updated.2019, depends on the names of genes in the data
g2m.genes <- union(cc.genes$g2m.genes, cc.genes.updated.2019$g2m.genes)
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

saveRDS(obj, "GSE152048_filt_all_seurat.RDS")

