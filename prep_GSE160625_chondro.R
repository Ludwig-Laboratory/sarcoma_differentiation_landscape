# Prepare data from GSE160625 into a single Seurat object
# GSE160625 citation: Wu et al. Nat Commun 2021
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160625

library(tidyverse)
library(Seurat)

#setwd("GSE160625_RAW") # specify directory with raw data downloaded from GEO


##### LOAD AND COMBINE #####

# get unique experiment sample names (of separate 10X datasets)
files = list.files()
data_names = files %>% str_extract(".+(?=_barcodes)") %>% na.omit %>% unique

# go through each unique sample to create a subdirectory copy of data
for (n in data_names) {
  # match files with same sample name
  files_i = files[str_detect(files,n)]
  new_files = str_extract(files_i, paste0("(?<=",n,"_).+")) # rename without prefix
  
  # create new dir, copy, and rename
  dir_i = n
  dir.create(dir_i)
  file.copy(files_i, dir_i)
  file.rename(file.path(dir_i,files_i), file.path(dir_i,new_files))
  file.rename(file.path(dir_i,"genes.tsv.gz"), file.path(dir_i,"features.tsv.gz")) # need to rename genes to features for Seurat
}

# load each sample with Seurat Read10X
obj_list = list()
for (i in 1:length(data_names)) {
  dir_i = data_names[i]
  data_i = Read10X(dir_i)
  obj_i = CreateSeuratObject(data_i, project=data_names[i])
  obj_list[[i]] = obj_i
}

# merge into one Seurat object
obj = merge(obj_list[[1]], obj_list[-1], add.cell.ids = data_names, project="GSE160625")

# save as .RDS
saveRDS(obj, "GSE160625_merged_seurat.RDS")


##### PROCESSING #####

# pre-processing and filtering
obj = readRDS("GSE160625_merged_seurat.RDS")

# add metadata (all we have is experimental condition/timepoint)
# for time, we set Cp to day 0, Scl = -6 days, iPSC = -12 days (per their chondrogenesis protocol)
obj$Condition = obj$orig.ident %>% str_extract("(?<=_).+")
obj$Condition = factor(obj$Condition, levels = obj$Condition %>% unique)
obj$Time = obj$Condition %>% str_extract("(?<=D).+") %>% as.numeric
obj$Time[obj$Condition == "Cp"] = 0
obj$Time[obj$Condition == "Scl"] = -6
obj$Time[obj$Condition == "hiPSC"] = -12


# perform QC filtering of low quality cells
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt", assay = 'RNA')

library(SingleCellExperiment)
library(scuttle)
library(scater)
obj.sce = as.SingleCellExperiment(obj)
discard.mito <- isOutlier(obj.sce$percent.mt, type="higher", batch=obj.sce$orig.ident)
discard.counts <- isOutlier(obj.sce$nCount_RNA, log = TRUE, type="lower", batch=obj.sce$orig.ident)
discard.features <- isOutlier(obj.sce$nFeature_RNA, log = TRUE, type="lower", batch=obj.sce$orig.ident)

# filter low quality cells based on: outliers with >3 mad in mito-mapped reads, <3 mad in nCount (total reads) or nFeature (unique genes mapped)
# also include a hard cutoff filter of mito<25%, nCount>1000, nFeature>500
obj.filt <- subset(obj, subset = !(discard.mito | discard.counts | discard.features) &
                       percent.mt < 25 & nCount_RNA > 1000 & nFeature_RNA > 500)

# some QC plots
plotColData(obj.sce, x="orig.ident", y="percent.mt", colour_by=I(discard.mito)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + labs(fill = 'Outlier') + ylab('% Mitochondrial Mapping')
plotColData(obj.sce, x="orig.ident", y="nCount_RNA", colour_by=I(discard.counts)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + ylab('Counts') + scale_y_log10()
plotColData(obj.sce, x="orig.ident", y="nFeature_RNA", colour_by=I(discard.features)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + ylab('Genes') + scale_y_log10()
RidgePlot(obj, group.by = "orig.ident", features = c("nCount_RNA")) + ggtitle("Counts")
RidgePlot(obj, group.by = "orig.ident", features = c("percent.mt")) + ggtitle("% Mitochondrial")
FeatureScatter(obj.filt, 'percent.mt', 'nFeature_RNA', pt.size = 0.1) + NoLegend() +
  xlab('% Mitochondrial Mapping') + ylab("Counts per cell") + theme(axis.text.y=element_text(size=6)) +
  xlim(0,100) + scale_y_log10(limits=c(500,10000))


##### SAVE DATA #####

# subset specific conditions within C59-treated differentiation timecourse
keep_cond = c("Cp", "C59_D7", "C59_D14", "C59_D28", "C59_D42")
obj.subset = obj.filt[,obj.filt$Condition %in% keep_cond]

# save all data and C59 subset data Seurat objects
saveRDS(obj.filt, "GSE160625_filtnorm_all_seurat.RDS")
saveRDS(obj.subset, "GSE160625_filtnorm_C59subset_seurat.RDS")

