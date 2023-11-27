# script for figure 2

library(Seurat)
library(harmony)
library(tidyverse)
library(ComplexHeatmap)


# Figure 2: Mesenchymal Transition Landscape

#Schematic of study design
#MTL UMAP: clusters, lineage, time (scaled to experiment)
#- No phase
#Marker genes: 2 per lineage



##### load data #####

# Load Seurat object from .rds file
# contains Osteo/Adipo + Chondro data, merged + SCTransform normalized by batch
datadir = "Data"
MTL_obj = readRDS( file.path(datadir, "MTL_OAC_SCT_qnorm", "MTL_OAC_merged_SCT.rds") )

# calculate scaled time as time divided by max for that experiment
max.times = MTL_obj@meta.data %>% select(Time,Lineage) %>%
  group_by(Lineage) %>% summarize(max.time = max(Time))
MTL_obj$time.scaled = MTL_obj$Time / max.times$max.time[match(MTL_obj$Lineage, max.times$Lineage)]



# run UMAP dimension reduction (PRE harmony batch correction)
MTL_obj = MTL_obj %>% ScaleData %>% RunPCA %>% RunUMAP(dims=1:30) # for umap visualization with DimPlot, FeaturePlot
DimPlot(MTL_obj, group.by="Condition", label=T, repel=T) + NoLegend()
DimPlot(MTL_obj, group.by="Lineage")
FeaturePlot(MTL_obj, "Time")

# find clusters
MTL_obj <- FindNeighbors(MTL_obj, dims = 1:10, reduction="harmony")
MTL_obj <- FindClusters(MTL_obj, resolution = 0.1)

# visualize before harmony
DimPlot(MTL_obj, group.by = 'seurat_clusters', reduction = 'umap')
DimPlot(MTL_obj, group.by = 'anno', reduction = 'umap', label=T) # previous labels from Osteo/Adipo data only
DimPlot(MTL_obj, group.by = 'Lineage', reduction = 'umap', label=T) + NoLegend()
DimPlot(MTL_obj, group.by = 'Phase', reduction = 'umap')
FeaturePlot(MTL_obj, features = "Time", reduction="umap")
FeaturePlot(MTL_obj, features = "time.scaled", reduction="umap")



##### harmony batch correction #####
# these parameters seem to work OK
#theta = 1, lambda = 0.5,  tau = 0
#theta = 1.5, lambda = 0.2,  tau = 200
#theta = 3, lambda = 0.5,  tau = 300 # selected for final figure
MTL_obj <- RunHarmony(MTL_obj, group.by.vars = 'Batch.ID', assay.use = 'SCT', dims.use = 1:50, 
                      max.iter.harmony = 10, reference_values = '2',
                      theta = 3, lambda = 0.5,  tau = 300,
                      #plot_convergence=T
                      ) %>%
  RunUMAP(dims = 1:30, reduction = 'harmony', reduction.name = 'umap_harmony')

# # loop over harmony parameters to optimize: 
# # save plots of the Harmony convergence and the resulting UMAP (two lists of plots for two pdfs)
# # make grid spanning theta, lambda over 5 values, 5 pages for tau
# theta_range = c(1,1.5,2,3,4)
# lambda_range = c(0.1,0.2,0.5,0.7,1)
# tau_range = c(0,50,100,200,300)
# umap_list1 = list()
# umap_list2 = list()
# ii = 1
# for (i_tau in 1:length(tau_range)) {
#   for (j_lambda in 1:length(lambda_range)) {
#     for (k_theta in 1:length(theta_range)) {
#       MTL_obj <- RunHarmony(MTL_obj, group.by.vars = 'Batch.ID', assay.use = 'SCT',
#                             max.iter.harmony = 40, dims.use = 1:50, reference_values = '2',
#                             theta = theta_range[k_theta],
#                             lambda = lambda_range[j_lambda],
#                             tau = tau_range[i_tau], verbose=F) %>%
#         RunUMAP(dims = 1:30, reduction = 'harmony', reduction.name = 'umap_harmony', verbose=F)
#       
#       print( sprintf("Plotting %i/125: tau=%i, lambda=%0.1f, theta=%0.1f",
#                      ii, tau_range[i_tau], lambda_range[j_lambda], theta_range[k_theta]) )
#       umap_list1[[ii]] = DimPlot(MTL_obj, group.by = 'Condition', reduction = 'umap_harmony',
#                                 label=T, repel=T) + NoLegend() +
#         ggtitle(bquote(list(lambda==.(lambda_range[j_lambda]),theta==.(theta_range[k_theta]))))
#       umap_list2[[ii]] = DimPlot(MTL_obj, group.by = 'Lineage', reduction = 'umap_harmony',label=T) + NoLegend() +
#         ggtitle(bquote(list(lambda==.(lambda_range[j_lambda]),theta==.(theta_range[k_theta]))))
#       ii=ii+1
#     }
#   }
# }
# library(ggpubr)
# pdf(file.path(resdir,"UMAP_harmony_range_10it_condition.pdf"), width=15, height=15)
# ggarrange(plotlist=umap_list1[1:25],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[1])))
# ggarrange(plotlist=umap_list1[26:50],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[2])))
# ggarrange(plotlist=umap_list1[51:75],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[3])))
# ggarrange(plotlist=umap_list1[76:100],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[4])))
# ggarrange(plotlist=umap_list1[101:125],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[5])))
# dev.off()
# pdf(file.path(resdir,"UMAP_harmony_range_10it_lineage.pdf"), width=15, height=15)
# ggarrange(plotlist=umap_list2[1:25],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[1])))
# ggarrange(plotlist=umap_list2[26:50],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[2])))
# ggarrange(plotlist=umap_list2[51:75],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[3])))
# ggarrange(plotlist=umap_list2[76:100],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[4])))
# ggarrange(plotlist=umap_list2[101:125],nrow=5,ncol=5) %>%
#   annotate_figure(top = text_grob(paste0("tau=",tau_range[5])))
# dev.off()



##### CLUSTERING (post harmony integration) #####
# find clusters after harmony correction
MTL_obj <- FindNeighbors(MTL_obj, k.param = 50, dims = 1:30, reduction="harmony", graph.name=c("harmony_nn","harmony_snn"))
MTL_obj <- FindClusters(MTL_obj, resolution = 0.3, graph.name="harmony_snn")
MTL_obj[["harmony_clusters"]] <- Idents(MTL_obj)

# add manual cluster label annotations
Idents(MTL_obj) <- MTL_obj[["harmony_clusters"]]
MTL_obj = RenameIdents(MTL_obj, `0` = "UD",
                       `1` = "A1",
                       `2` = "O2",
                       `3` = "A3",
                       `4` = "MSC-H",
                       `5` = "O1",
                       `6` = "C3",
                       `7` = "C2",
                       `8` = "C1",
                       `9` = "CP",
                       `10` = "MSC-L",
                       `11` = "A2",
                       `12` = "A4")
levels(MTL_obj) <- c("UD", "MSC-H", "MSC-L",
                     "O1", "O2", "A1", "A2", "A3", "A4", "CP", "C1", "C2", "C3")
MTL_obj$cluster_annotation = Idents(MTL_obj)


# UMAPs: clusters, lineage, time (scaled)
DimPlot(MTL_obj, reduction = 'umap_harmony', label=T) + theme(legend.position = "bottom", aspect.ratio = 1,
                                                              axis.title.x = element_blank(),
                                                              axis.title.y = element_blank())
DimPlot(MTL_obj, group.by = 'Lineage', reduction = 'umap_harmony', label=T) + NoLegend() + theme(aspect.ratio = 1,
                                                                                                 axis.title.x = element_blank(),
                                                                                                 axis.title.y = element_blank())
FeaturePlot(MTL_obj, features = "time.scaled", reduction="umap_harmony") +
  ggtitle("Time (scaled)") +
  theme(aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) & viridis::scale_color_viridis()
#FeaturePlot(MTL_obj, features = "Time", reduction="umap_harmony") & viridis::scale_color_viridis()
#DimPlot(MTL_obj, group.by = 'anno', reduction = 'umap_harmony', label=T)
#DimPlot(MTL_obj, group.by = 'Phase', reduction = 'umap_harmony')


# counts of lineage in each cluster
MTL_obj@meta.data %>% select(Lineage, cluster_annotation) %>%
  group_by(cluster_annotation) %>% summarize(adipo = sum(Lineage=="Adipo"),#/n(),
                                           chondro = sum(Lineage=="Chondro"),
                                           osteo = sum(Lineage=="Osteo")) %>%
  reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=cluster_annotation, y=value, fill=variable), stat="identity") +
  ylab("count")

# fraction of timepoints in each cluster
MTL_obj@meta.data %>% select(Time, cluster_annotation) %>%
  group_by(cluster_annotation, Time) %>% summarize(n=n()) %>%
  filter(cluster_annotation == "A3")




##### MARKER GENE ANALYSIS #####
# FindAllMarkers and save
# identify differentially expressed gene markers
MTL_obj = PrepSCTFindMarkers(MTL_obj) # may take a while
markers.all = FindAllMarkers(MTL_obj)
markers.all %>% group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>% 
  mutate(clust_genes = paste0(gene, collapse = ";")) %>%
  select(cluster, clust_genes) %>% distinct

resdir = "Results"

# write markers to a csv
write_csv(markers.all, file.path(resdir, "MTL_cluster_markers.csv"))

markers.pos = FindAllMarkers(MTL_obj, only.pos = T)
write_csv(markers.pos, file.path(resdir, "MTL_cluster_markers_pos.csv"))



# Marker Gene plots:
# violins, 2 genes per lineage

# violin plots (using clusters)
genes = c("COL1A1", "COL8A1","MGP","DCN", # osteo markers: also RUNX2, LUM, SPARC, POSTN
          "FABP4", "ADIPOQ", "ACACB", "PLIN1", # adipo markers: also PPARG, CEBPA, LPL, PLIN4
          "COL2A1", "SOX9", "ACAN", "PTHLH") # chondro markers:
#VlnPlot(MTL_obj, features=genes, group.by="harmony_clusters", ncol=4)
VlnPlot(MTL_obj, features=genes, group.by="harmony_clusters", ncol=4, pt.size=0)


# osteo markers
genes = c("COL1A1", "COL8A1","MGP","DCN", "RUNX2", "LUM", "SPARC", "POSTN")
genes = c("COL1A1", "COL8A1","MGP","DCN", "RUNX2", "LUM", "SPARC", "POSTN")
VlnPlot(MTL_obj, features=genes, group.by="harmony_clusters", ncol=4, pt.size=0)

# adipo markers
genes = c("FABP4", "ADIPOQ", "ACACB", "PLIN1", "PPARG", "CEBPA", "LPL", "PLIN4")
VlnPlot(MTL_obj, features=genes, group.by="harmony_clusters", ncol=4, pt.size=0)

# chondro markers
genes = c("COL2A1", "SOX9", "ACAN", "PTHLH")
VlnPlot(MTL_obj, features=genes, group.by="harmony_clusters", ncol=2, pt.size=0)


# specified markers
genes = c("CCND1", "FGF2", "THY1", "NT5E", "TUBA1C", "VIM", "PRDX1", "HAS1",
          "C1R", "DCN", "FGF7", "CDH11", "PDGFRA", "TWIST1", "CD44", "COL6A3",
          "RUNX2", "TGFBI", "WNT5A", "HAS2", "PPARG", "FABP4", "ADIPOQ", "CEBPA",
          "ACTA2", "CRYAB", "COL1A1", "ADAMTS7", "CXADR", "IGDCC3", "MDK", "FGF12",
          "HEY1", "ASPN", "GAS2", "FOS", "SOX9", "PTH1R", "COL2A1", "ACAN",
          "KLF4", "IFITM1", "EIF1", "CYTL1", "ELN", "POSTN", "NOX4", "FBLN5")
VlnPlot(MTL_obj, features=genes, group.by="harmony_clusters", ncol=8, pt.size=0)

# subselected markers
genes = c("LEPR","CLIC3","TAGLN",
          "IGFBP6", "FOXC2", "GLIPR1", "CCND1",
          "SAA1", "SAA2", "FBLN1", "CHI3L1", "THBS2") # osteo markers
genes = c("MT1X", "MT1E", "MT1M", "MT2A",
          "WNT5A", "PAPPA", "FTH1",
          "LUM", "COMP", "COL3A1", "COL6A1", "DCN") # adipo markers
genes = c("CLIC3", "FABP4", "COL2A1", "IGFBP5")
FeaturePlot(MTL_obj, genes, reduction="umap_harmony")

# final selected markers
genes = c("CLIC3","LEPR", #osteo
          "MT1X","FABP4", #adipo
          "MDK", "COL2A1") #chondro
VlnPlot(MTL_obj, features=genes, ncol=2, pt.size=0) & labs(x=NULL)



# save merged data with harmony coordinates
saveRDS(MTL_obj, file.path(datadir,"MTL_OAC_SCT_qnorm/MTL_OAC_harmony.rds"))
MTL_obj = readRDS(file.path(datadir,"MTL_OAC_SCT_qnorm/MTL_OAC_harmony.rds"))



#### SUPPLEMENTAL FIGURE 1: YAP/TAZ in MTL

# markers in MSC-H vs MSC-L
markers.msc = FindMarkers(MTL_obj, ident.1 = "MSC-H", ident.2 = "MSC-L")
yaptaz = c("TEAD1", "TEAD2", "CTGF", "CYR61", "IGFBP5", "ANKRD1")

# demonstrate YAP1/TAZ transcriptional stimulation
# genes from MSigDB C2: Curated Canonical Pathways (CP:Reactome): REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION
# https://www.gsea-msigdb.org/gsea/msigdb/cards/REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION
# and another reference: CTGF (alias CCN2), CYR61, IGFBP5, ANKRD1 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5190017/)
yaptaz = c("YAP1", "TAZ", "WWTR1",
           "CTGF", "GATA4", "HIPK1", "HIPK2", "KAT2B", "NKX2-5",
           "TBX5", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
           "CYR61", "IGFBP5", "ANKRD1") %>% sort

# score pathway in MTL
MTL_obj <- AddModuleScore(MTL_obj, list(yaptaz), name = "YAPTAZ")

# plot UMAP colored by YAP/TAZ and Violins of each cluster
p1 <- FeaturePlot(MTL_obj, "YAPTAZ1", reduction="umap_harmony", label=T) +
  theme(aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + ggtitle("YAP/TAZ signature")
p2 <- VlnPlot(MTL_obj, "YAPTAZ1", pt.size = 0,
              fill.by = "ident") +
  theme(axis.title.x = element_blank(), legend.position = "none") + ggtitle("YAP/TAZ signature")
ggpubr::ggarrange(p1,p2,widths = c(1.2, 1), align = "h")

# dotplot of each gene
p3 <- DotPlot(MTL_obj, features = yaptaz) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y = element_blank()) +
  labs(x = "YAP/TAZ Genes") +
  scale_y_discrete(limits=rev)
p3
