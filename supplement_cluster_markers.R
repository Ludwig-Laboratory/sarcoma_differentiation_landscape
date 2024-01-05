# supplemental figure of MTL cluster analysis
# we will examine marker genes (already determined with Seurat's FindAllMarkers)
# 1. Heatmap of selected marker genes (to illustrate cluster identities with known molecular pathways)
# 2. Pathway enrichment plots of marker genes in each cluster

library(tidyverse)
library(Seurat)

# specify directory to save enrichment results
resdir = "Results/enrichment"

##### LOAD DATA #####
# load gene expression data from MTL (osteo+adipo+chondro, SCT normalized)
# already with annotated clusters
MTL_obj = readRDS( "Data/MTL_OAC_SCT_qnorm/MTL_OAC_harmony.rds" )

# load in cluster marker genes, pre-computed in script figure2.R
cluster_markers = read_csv( 'Results/MTL_cluster_markers_pos.csv' )


##### PANEL A: MARKER GENE HEATMAP #####

# heatmap of selected marker genes
selected_genes = c("ENG", "NT5E", "PRRX1", "THY1", # MSC markers
                   "YAP1", "WWTR1", "ANKRD1", "CTGF", "CYR61", "IGFBP5", "TEAD1", # YAP/TAZ
                   "COL3A1", "COL6A1", "COMP", "DCN", "LUM", # extracellular matrix
                   "COL1A1", "COL1A2", "COL8A1", "ELN", # ECM again
                   "WNT5A", "PAPPA", "FTH1", # adipogenesis
                   "MT1X", "MT1E", "MT1M", "MT2A", # metallothioneins
                   "ACACB", "ADIPOQ", "APOE", "FABP4", "G0S2", "FABP5", "LPL", "PLIN4", "PLIN1", # adipocyte
                   "WISP2", "AXL", "DKK1", "HHIP", # early osteogenesis
                   "CCND1", "FOXC2", "GLIPR1", # GLI targets
                   "FN1", "GSN", "SERPINE2", # osteoblast secretome
                   "CHI3L1", "FBLN1", "SAA1", "THBS2", # osteoblast secreteome
                   "FOS", "CTNNB1", # mechanosensing
                   "IGFBP2", "IGFBP6", "IGFBP4", "IGFBP7", "IGF2", # IGFPs
                   "SOX2", "SOX4", "SOX6", # early chondrogenic potential
                   "COL9A1", "MATN4", "SOX9", # chondrogenic
                   "ACAN", "COL2A1", # later chondro
                   "EPYC", "FRZB", "LECT1" # frizzled-related
)
# annotate gene classes for each gene
gene_class = c(rep("MSC", 4),
               rep("YAP/TAZ", 7),
               rep("Extracellular Matrix", 9),
               rep("Adipogenesis", 3),
               rep("Metallothioneins", 4),
               rep("Adipocyte", 9),
               rep("Early osteogenesis", 4),
               rep("GLI targets", 3),
               rep("Osteoblast secretome", 7),
               rep("Mechanosensing", 2),
               rep("IGFPs", 5),
               rep("Early Chondrogenesis", 3),
               rep("Later Chondrogenesis", 5),
               rep("Frizzled-related", 3))
gene_class = factor(gene_class, levels=unique(gene_class))

expr = MTL_obj@assays$SCT@data # using SCTransform normalized expression
expr_subset = expr[selected_genes,] %>% as.matrix
expr_scaled = expr_subset/apply(expr_subset,1,max)

col_df = data.frame(anno = MTL_obj$cluster_annotation) %>%
  mutate(lineage = if_else(anno %in% c("UD","MSC-H","MSC-L","CP"), "Undifferentiated",
                           if_else(anno %in% c("O1","O2"), "Osteogenic",
                                   if_else(anno %in% c("A1","A2","A3","A4"), "Adipogenic",
                                           if_else(anno %in% c("C1","C2","C3"), "Chondrogenic", "unknown")))))
ta = HeatmapAnnotation(which = "column",
                       lineage = col_df$lineage,
                       col = list(lineage = c("Adipogenic" = "#F8766D", "Chondrogenic" = "#00BA37", "Osteogenic" = "#619CFF", "Undifferentiated" = "#808080")))

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("#009FFF", "#FFFF00", "#FF5100"))

pdf( file.path(resdir, "MTL_ClusterMarkers_Heatmap.pdf"), width=20, height=10)
Heatmap(expr_scaled, name="scaled norm. expr.", col=col_fun,
        row_split = gene_class, cluster_row_slices = F, # cluster_rows = F,
        row_title_rot = 30, row_title_side = "right", show_row_dend = F,
        top_annotation = ta,
        column_split = MTL_obj$cluster_annotation %>% fct_recode(MSCH="MSC-H",MSCL="MSC-L"),
        cluster_column_slices = F,
        show_column_dend = F, show_column_names =F)
dev.off()



##### PANEL B: PATHWAY ENRICHMENT #####

# GO enrichment based on list of marker genes of each cluster
# use clusterProfiler enricher
library(clusterProfiler)
library(enrichplot)


# first convert HGNC gene symbol to Entrez gene ID
library(org.Hs.eg.db)
genes_HGNC = cluster_markers$gene # cluster marker genes as HGNC symbols
genes_entrez = mapIds(org.Hs.eg.db,
                      keys = genes_HGNC,
                      keytype = "SYMBOL",
                      column = "ENTREZID",
                      multiVals = "first")


# compute GO enrichment of each cluster
clusters = unique(cluster_markers$cluster); nclust = length(clusters)

ego_list = list() # store each cluster enrichment results in list
plot_list = list() # create enrichment dotplots of each cluster
for (i in 1:nclust) {
  cind = which(cluster_markers$cluster == clusters[i])
  cluster_genes = genes_entrez[cind] # entrezid of marker genes to enrich
  
  print(paste("Cluster",i,": ngenes =",length(cluster_genes)))
  ego = enrichGO( gene = cluster_genes,
                  universe = genes_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = T)
  ego_list[[i]] = ego
}
for (i in 1:nclust) {
  plot_list[[i]] = dotplot(ego_list[[i]], font.size = 8, showCategory=10) +
    ggtitle(sprintf( "Cluster %i:%s (n=%i genes)", i, clusters[i], sum(cluster_markers$cluster==clusters[i]) )) +
    scale_color_gradient(low="red", high="blue", limits = c(NA, 0.01),
                         trans =  "log10",
                         labels = function(x) sprintf("%.1e",x)) +
    theme(axis.text.y = element_text(size=6),
          plot.title = element_text(size = 10))
}
sapply(ego_list,dim)

# combine all clusters into one plot
p_enrich = ggpubr::ggarrange(plotlist=c(plot_list[1:9],list(NULL),plot_list[10:13]), nrow=3, ncol=5,
                             common.legend=T, legend="right")
pdf( file.path(resdir, "MTL_Cluster_EnrichGO.pdf"), width=20, height=10)
p_enrich
dev.off()

# save each cluster enrichment to a excel worksheet page
library(xlsx)
library(rJava)
options(java.parameters = "-Xmx1000m")
filename = "Results/enrichment/MTLcluster_enrichGO.xlsx"
# first page is the cluster marker genes table
write.xlsx(as.data.frame(cluster_markers),
           file=filename, sheetName = "cluster_marker_genes", row.names=F)
# each subsequent page is enrichment results of each cluster
for (i in 1:nclust) {
  invisible(gc()); .jcall("java/lang/System", method = "gc") # this line clears memory between each page (or you might get a heap space out of memory error)
  ego_list[[i]] %>% as.data.frame %>%
    arrange(p.adjust) %>%
    write.xlsx(file=filename,
               sheetName=paste0("cluster",i,"_",clusters[i],"_enrichGO"), row.names=FALSE,
               append=T) # append pages to xlsx sheet
}


# # CellMarker
# # from : http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx'
# cell_marker = readxl::read_xlsx( file.path(resdir,"Cell_marker_Human.xlsx"))
# 
# # # with Entrez ID
# # cells <- cell_marker %>%
# #   dplyr::select(cell_name, GeneID) %>%
# #   dplyr::mutate(GeneID = strsplit(GeneID, ', ')) %>%
# #   tidyr::unnest()
# # with HGNC gene symbol
# cells <- cell_marker %>%
#   dplyr::select(cell_name, Symbol) %>%
#   dplyr::mutate(Symbol = strsplit(Symbol, ', ')) %>%
#   tidyr::unnest()
# 
# cell_marker %>% select(tissue_class, cell_name) %>% # tissue_class tissue_name cell_type cell_name
#   distinct() %>% dim
# 
# ego_list2 = list()
# for (i in 1:nclust) {
#   cind = which(cluster_markers$cluster == clusters[i])
#   cluster_genes = genes_HGNC[cind] # entrezid of marker genes to enrich
#   
#   print(paste("Cluster",i,": ngenes =",length(cluster_genes)))
#   ego = enricher(cluster_genes, TERM2GENE = cells)
#   ego_list2[[i]] = ego
# }
# 
# filename = file.path(resdir,"cluster_clustermarker.xlsx")
# write.xlsx(cluster_markers,
#            file=filename, sheetName = "cluster_marker_genes", row.names=F)
# for (i in 1:nclust) {
#   ego_list2[[i]] %>% as.data.frame %>%
#     write.xlsx(file=filename,
#                sheetName=paste0("cluster",i,"_",clusters[i],"_enrichGO"), row.names=FALSE,
#                append=T) # append pages to xlsx sheet
# }
