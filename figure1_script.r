#Load Libraries
library(Seurat)
library(ggplotify)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)

#Make figure 1 for the MTL
#Panel A and B: UMAP of Osteosarcoma (OS) PDX single-cell data, A is sample ID and B is cluster
#Panel C: Heatmap of gene signatures for OS PDX
#Panel D: Heatmap of differentially expressed genes
#Panel E: Dot plot of pathway analysis

#Load the OS PDX data set
OS_PDX_Warm <- readRDS( 'OS_PDX_Warm.rds')

##Generating UMAPs for MTL Paper
require(ggrepel)
#make dummy umap to take the labels
p <- DimPlot(OS_PDX_Warm, group.by = 'seurat_clusters', pt.size = NA, label = TRUE, repel = TRUE) + NoLegend() 
p <- ggplot_build(p)
df <- p[["data"]][[2]]
rm(p)

umap_plots <- list() #initialize a list of umap

#make a plot with the patients identified
umap_plots[[1]] <- DimPlot(OS_PDX_Warm, reduction = "umap", group.by = "lab_id", label = T, repel = T) +  guides(color=guide_legend(title="Lab ID", nrow = 1, override.aes = list(size = 3)))

#make a plot with cells identified from clusters
umap_plots[[2]] <- DimPlot(OS_PDX_Warm, reduction = "umap", group.by = "seurat_clusters", label = FALSE, label.size = 4, repel = TRUE) + geom_label_repel(df, mapping = aes(x=x,y=y, label = label, fontface = 'bold'), fill = alpha(c("white"),3/5), label.size= NA , box.padding = 0, label.r = 0.5) + guides(color=guide_legend(title="Cluster", nrow = 2, override.aes = list(size = 2))) 

umap_plots <-
  lapply(umap_plots, function(x) {
    x + theme_classic() + theme(
      legend.position = 'top',
      legend.justification = 'left',
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = 'bold'),
      legend.spacing.x= unit(0.04, 'mm'),
      axis.text = element_blank(),
      plot.title = element_text(size = 10)
    ) + ggtitle('')
  })

#View UMAPs
umap_plots


#Heatmap from Signature Genes
cell_sigs <- list(Osteoblastic = c('RUNX2', 'SPP1', 'COL1A1', 'CDH11', 
                                   'IBSP','PTH1R'),
                  Chondroblastic = c('COL2A1', 'SOX9', 'ACAN','PTHLH'),
                  Fibroblastic = c('POSTN', 'PDGFRA','ACTA2',
                                   'THBS1'),
                  Adipocytic = c('PLIN2', 'FABP4', 'LPL'),
                  EMT  = c('SNAI1', 'SNAI2', 'TWIST1', 'VIM'),
                  Stemness = c('THY1', 'NT5E','MYC','NES', 'KLF4'),
                  'Cycling Cells' = c('MKI67', 'TOP2A', 'CENPF'))

downsample <- subset(OS_PDX_Warm, downsample = 1000) #speeds up the heatmap
metadata <- downsample@meta.data

require(scales)
mat_scaled <- downsample@assays$RNA@data[unlist(cell_sigs),]
mat_scaled <- apply(mat_scaled, 1, function(x) rescale(x, c(0,1))) %>% t()

gene_anno <- data.frame(genes = unlist(cell_sigs)) #collapses the list into a dataframe

gene_anno$signature <- ''
for (x in names(cell_sigs)){
  gene_anno[grep(pattern = paste0('^', x), x = rownames(gene_anno)), 'signature'] <- x
}
rownames(gene_anno) <- gene_anno$genes

#Generate the colors for each group
group_colors <- setNames(brewer.pal(n =length(
  unique(names(cell_sigs))), name = "Set3"), unique(names(cell_sigs)))

cluster_colors <- setNames(gg_color_hue(n =length(
  unique(metadata$lab_id))), unique(metadata$lab_id))

annotation_colors <- list(signature = group_colors,
                          lab_id = cluster_colors)

#Annotate each gene based on which group they belong in

ra_1 = rowAnnotation(
  df =  gene_anno[rownames(mat_scaled), 'signature', drop = FALSE],
  col = annotation_colors,
  annotation_label = '',
  border = TRUE,
  show_legend = F,
  simple_anno_size  = unit(0.25, 'cm'),
  annotation_legend_param = list(signature = list(
    nrow = 1, title_position = "lefttop"
  ))
  
) 

lgd_anno = Legend(labels = names(cluster_colors), 
                  title = "PDX", ncol = 2,
                  title_position = 'topleft',
                  direction = c( "horizontal"),
                  legend_gp = gpar(fill = cluster_colors))

#annotate each cell based on which sample they are from
ha_1 = HeatmapAnnotation(
  df =  metadata[, 'lab_id', drop = FALSE],
  col = annotation_colors,
  annotation_label = 'Lab ID',
  border = F,
  show_legend = F,
  simple_anno_size  = unit(0.25, 'cm'),
  annotation_legend_param = list(lab_id = list(
    nrow = 1, title_position = "lefttop"
  ))
)

sig_hm <- Heatmap(
  mat_scaled,
  col = f1,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  row_split = gene_anno[rownames(mat_scaled), 'signature', drop = FALSE], #split based on signatures
  cluster_row_slices = F,
  show_row_names = TRUE,
  show_column_names = FALSE,
  right_annotation = ra_1,
  top_annotation = ha_1,
  heatmap_height = unit(6, 'in'),
  heatmap_width = unit(6, 'in'),
  column_names_rot = 90,
  column_title_rot = 0,
  column_title_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 10),
  row_title_side = 'right', row_title_rot = 0,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Normalized Expression', 
                              direction = "horizontal", legend_width = unit(3, 'cm')))

#Save Heatmap 
sig_hm_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(
  sig_hm,
  merge_legend = T ,
  annotation_legend_list = list(lgd_anno),
  heatmap_legend_side = c("bottom"),
  annotation_legend_side = c("bottom")
))

#Find DEGs

#plan("multiprocess", workers = 8)
#OS_PDX_Warm_markers <- FindAllMarkers(OS_PDX_Warm, test.use = 'LR', only.pos = T, return.thresh = 0.05)
#plan("sequential")
#saveRDS(OS_PDX_Warm_markers, 'OS_PDX_Warm markers.rds')

#Load the DEGS
OS_PDX_Warm_markers <- readRDS('OS_PDX_Warm markers.rds')

#Prepare DEG Heatmap
#use the top 3 markers from each cluster
top_markers <- OS_PDX_Warm_markers %>% group_by(cluster) %>% top_n(3, avg_log2FC) %>% pull(gene) %>% unique() 
OS_PDX_Warm@active.ident <- OS_PDX_Warm$seurat_clusters
OS_PDX_Warm_avg <- AverageExpression(OS_PDX_Warm, features = top_markers)

metadata <- data.frame(clusters = colnames(OS_PDX_Warm_avg[["RNA"]]))
rownames(metadata) <- metadata$clusters

require(scales) #re-scale the data
mat_scaled <- OS_PDX_Warm_avg[["RNA"]][top_markers,]
mat_scaled <- apply(mat_scaled, 1, function(x) rescale(x, c(0,1))) %>% t()


cluster_colors <- setNames(gg_color_hue(n =length(
  unique(metadata$cluster))), unique(metadata$cluster)) #generating the colors for each cluster

annotation_color <- list(
  clusters = cluster_colors)

# Make Stacked Bar Plot that goves over heatmap
require(tidyr)

table.cluster.experiment <- table(OS_PDX_Warm$lab_id, OS_PDX_Warm$seurat_clusters) #create table with lab_id and clusters

freq.cluster.experiment <- table.cluster.experiment / rowSums(table.cluster.experiment) 
freq.experiment.cluster <- sweep(table.cluster.experiment, 2, colSums(table.cluster.experiment), '/') 
freq_df <- as.data.frame.matrix(freq.experiment.cluster) %>% t() *100

lab_id_colors <- setNames(gg_color_hue(n =length(
  colnames(freq_df))), colnames(freq_df))

#Create DEG Heatmap
require(ComplexHeatmap)
f1 = colorRamp2(c(0,1), c("white", "darkblue"), transparency = 0, space = "LAB") #create color scale

a_1 = HeatmapAnnotation(
  df =  metadata[, 'clusters', drop = FALSE],
  col = annotation_color,
  annotation_label = 'Cluster',
  border = TRUE,
  show_legend = F,
  simple_anno_size  = unit(0.25, 'cm'),
  annotation_legend_param = list(clusters = list(
    nrow = 1, title_position = "lefttop"
  ))
)

cluster_colors_ = right_join(data.frame(names = metadata$cluster), data.frame(cluster_colors, names = names(cluster_colors)))

ha_2 = HeatmapAnnotation(
  rn = anno_text(metadata$cluster,
                 rot = 0,
                 gp = gpar(fontsize = 6, border = "black", fill = cluster_colors_$cluster_colors),
                 height = unit(4, 'mm'),
                 just = 'center',
                 location = unit(0.5, 'npc')
  ))

#Stacked Bar plot of lab_id frequencies in the clusters
ha_1 = HeatmapAnnotation('Pct. OS PDX' = anno_barplot(freq_df, bar_width = 1, 
                                                      gp = gpar(fill = lab_id_colors) ),
                         show_annotation_name = TRUE, 
                         annotation_name_side = 'left', height = unit(15, 'mm'),
                         annotation_name_gp = gpar(fontsize = 10)) 

lgd_anno = Legend(labels = names(lab_id_colors), 
                  title = "PDX", ncol = 2,
                  title_position = 'topleft',
                  direction = c( "horizontal"),
                  legend_gp = gpar(fill = lab_id_colors))


deg_hm <- Heatmap(
  mat_scaled,
  col = f1,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = fill))}, #adds gray outlines
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  cluster_row_slices = F,
  cluster_column_slices = F,
  top_annotation =  c(ha_1, ha_2),
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_height = unit(5.5, 'in'),
  heatmap_width = unit(4, 'in'),
  column_names_rot = 90,
  column_title_rot = 0,
  column_title_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 10),
  row_title_side = 'right', row_title_rot = 0,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Normalized Expression', 
                              direction = "horizontal", legend_width = unit(3, 'cm'))
)

#Convert DEG Heatmap to Grob
deg_hm_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(
  deg_hm,
  merge_legend = T ,
  annotation_legend_list = list(lgd_anno),
  heatmap_legend_side = c("bottom"),
))


#Pathway Analysis using clusterProfiler
require(clusterProfiler)
require(msigdbr)
require(org.Hs.eg.db)

#download the msigdb hallmark data set
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
head(h_gene_sets)

#Select top 100 DEGs for each cluster
top_markers <- OS_PDX_Warm_markers %>% group_by(cluster) %>% top_n(100, avg_log2FC)

#convert from gene symbol to ENTREZID
top_markers$ENTREZID <- AnnotationDbi::select(org.Hs.eg.db, 
                                              keys = top_markers$gene,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "SYMBOL")$ENTREZID

compKEGG <- compareCluster(geneCluster   =  ENTREZID~cluster,
                           data = top_markers,
                           fun = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
)

pathways_dot_plot <- dotplot(compKEGG, showCategory = 1, font.size = 13,) 

pathways_dot_plot


#Creating the final figure 1
top_left <- plot_grid(plotlist = umap_plots, nrow = 2, labels =  c('A', 'B')) #UMAP and signature heatmap
top <- plot_grid(top_left, sig_hm_grob, ncol = 2, labels = c('', 'C')) 

bot <- plot_grid(deg_hm_grob, pathways_dot_plot, 
                 rel_widths = c(1,1.6), labels = c('D', 'E')) #added deg heatmap and pathway analysis plot

tiff('OS_PDX_MTL_Figure_1.tiff' , width = 12, height = 15, res= 150, units = 'in')
plot_grid(top,bot, nrow=2, rel_heights = c(1,0.8))
dev.off()

svg('OS_PDX_MTL_Figure_1.svg' , width = 12, height = 15)
plot_grid(top,bot, nrow=2, rel_heights = c(1,0.8))
dev.off()

