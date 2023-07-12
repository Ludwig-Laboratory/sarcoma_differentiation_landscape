# script to produce figures 3 and 4

library(tidyverse)
library(ComplexHeatmap)


# load in NMF data (output from MATLAB analysis script)

# Gene-Archetype correlations, Composition % and Order for PDX and OS
datadir = "Results/MTL_Archetype_Results"
gene_cor = read_csv( file.path(datadir, "T_Gene_Ordered.csv") )
comp_pdx = read_csv( file.path(datadir, "T_Composition_PDX.csv") )
comp_os = read_csv( file.path(datadir, "T_Composition_OS.csv") )
order_pdx = read_csv( file.path(datadir, "T_PDX_Order.csv") )
order_os = read_csv( file.path(datadir, "T_OS_Order.csv") )

# Archetype scores for PDX and OS
arch_mtl = read_csv( file.path(datadir, "MTL_Arch.csv"), col_names = F)
arch_pdx = read_csv( file.path(datadir, "PDX_Arch.csv"), col_names = F)
arch_os = read_csv( file.path(datadir, "OS_Arch.csv"), col_names = F )
mtl_condition = read_csv( file.path(datadir, "MTL_Cell_Label.csv"), col_names = F)


# cell metadata for each dataset
meta_mtl = read_delim("Data/MTL_OAC_SCT_qnorm/cell_metadata.txt")
meta_pdx = read_delim("Data/OSPDX_SCT_qnorm/cell_metadata.txt")
meta_os = read_delim("Data/OS11_SCT_qnorm/cell_metadata.txt")

# rename MDA-SA98-TIS02 to SA98
meta_pdx$lab_id[meta_pdx$lab_id == "MDA-SA98-TIS02"] = "SA98"

# set numeric order of OS11 tumor IDs
meta_os = meta_os %>%
  mutate(tumor_id = factor(tumor_id,
                           levels = c("BC2", "BC3", "BC5", "BC6", "BC10", "BC11",
                                      "BC16", "BC17", "BC20", "BC21", "BC22")))


# check cell barcode order
# PDX order matches original data
all(meta_pdx$barcode == order_pdx$Barcode)

# OS does not match
all(meta_os$barcode == order_os$Barcode)

# reorder OS data
matchind = match(meta_os$barcode, order_os$Barcode)
arch_osm = arch_os[,matchind]


# filter out OS cell subpopulations based on clustering
os_clusters = read_csv("Data/OS11_SCT_qnorm/OS11_barcode_clusters.csv") %>%
  mutate(cluster_annotation = factor(cluster_annotation,
                                     levels=c("OB", "CB", "MSC", "Prolif", "OC", "Myel", "TIL", "Endo")))
all(os_clusters$barcode == meta_os$barcode)

# lastly, load PDX cluster labels
pdx_clusters = read_csv("Results/PDX_clusters/OS_PDX_clusters.csv") %>%
  mutate(seurat_clusters = factor(seurat_clusters))
all(pdx_clusters$barcode == meta_pdx$barcode)


# output directory for plots
outdir = "Results/figures/subfigures"



# figure 3A: MTL archetype timecourses
mtl_line = mtl_condition$X1 %>% str_sub(1,1)
mtl_time = mtl_condition$X1 %>% str_sub(2)
mtl_time[mtl_time == "p5"] = 0.5
mtl_time = as.numeric(mtl_time)

plot_list = list()
for (i in 1:12) {
  if (i %in% c(1:7)) {lim = 21} else {lim=42}
  plot_list[[i]] = data.frame(arch= as.numeric(arch_mtl[i,]), l = mtl_line, t = mtl_time) %>%
    group_by(l,t) %>%
    summarize(m = mean(arch), s = sd(arch), n=n()) %>%
    mutate(sem = s/sqrt(n)) %>%
    ggplot() + geom_line(aes(x=t, y=m, color=l)) +
    geom_errorbar(aes(x=t, ymin=m-sem, ymax=m+sem, color=l)) +
    theme_classic() + theme(legend.position = "none", aspect.ratio = 1) +
    labs(x="Time", y=paste0("Archetype ",i)) + coord_cartesian(xlim=c(0, lim))
}
pdf( file.path(outdir, "MTL_Archetype_Timecourses.pdf"),
               width = 12, height=4)
ggpubr::ggarrange(plotlist = plot_list, ncol=6, nrow=2)
dev.off()

# figure 3B: heatmap of gene archetype correlations
x = as.matrix(gene_cor[,-1]) %>% t
colnames(x) = gene_cor$Gene
pdf( file.path(outdir, "MTL_Archetype_Gene_Correlations.pdf"),
     width=12, height=4)
Heatmap(x, name = "Correlation",
        cluster_rows = F, row_names_side = "left",
        row_title = "Archetype", row_title_side = "left",
        cluster_columns = F,
        column_title = "Gene", column_title_side = "bottom")
dev.off()



# heatmap 4A: PDX archetype scores
x = as.matrix(arch_pdx)
rownames(x) = 1:12
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#009FFF", "#FFFF00", "#FF5100"))
gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n]}
cluster_colors = gg_color_hue(nlevels(pdx_clusters$seurat_clusters)) %>% setNames(levels(pdx_clusters$seurat_clusters))
tumor_colors = gg_color_hue(3) %>% setNames(unique(meta_pdx$lab_id))

pdf( file.path(outdir, "PDX_Arch_Heatmap.pdf"),
     width=10, height=4)
Heatmap(x, name = "NMF Score", col = col_fun,
        cluster_rows = F, row_names_side = "left",
        row_title = "Archetype", row_title_side = "left",
        column_split = meta_pdx$lab_id, cluster_column_slices = F,
        column_title = "%s", show_column_names = F, show_column_dend = F,
        #top_annotation = HeatmapAnnotation(PDX = meta_pdx$lab_id,
        #                                   cluster = pdx_clusters$seurat_clusters,
        #                                   col = list(cluster=cluster_colors,
        #                                              PDX=tumor_colors))
        )
dev.off()



# figure 4B: PDX composition plot

# recompute compositions
comp_pdx2 = data.frame(Archetype = apply(arch_pdx, 2, which.max), id = meta_pdx$lab_id) %>%
  mutate(Archetype = factor(Archetype, levels = c(1:12))) %>%
  group_by(id, Archetype) %>% summarize(n = n()) %>%
  reshape2::dcast(Archetype ~ id)
comp_pdx2[is.na(comp_pdx2)] = 0
comp_pdx2[,-1] = lapply(comp_pdx2[,-1], function(x) x/sum(x))


# save composition bar plot with grouped "Other" archetypes: 1,4,7,8,11,12
pdf( file.path(outdir, "PDX_Arch_Compositions.pdf"),
     width=4, height=4)
# sum "Other" archetypes: 1,4,7,8,11,12
comp_pdx2 %>% mutate(other = Archetype %in% c(1,4,7,8,11,12)) %>%
  group_by(other) %>% summarise_if(is.numeric,sum) %>% filter(other==T) %>%
  rename(Archetype = other) %>% mutate(Archetype = "Other") %>%
  rbind(comp_pdx2 %>% filter( !(Archetype %in% c(1,4,7,8,11,12))),.) %>%
  reshape2::melt() %>%
  ggplot() + geom_col(aes(y=variable, x=value, fill = Archetype),
                      position = position_fill(reverse = TRUE)) +
  scale_y_discrete(limits=rev) +
  theme_classic() + ylab("Tumor ID") + xlab("Fraction cells")
dev.off()



# heatmap 4C: OS11 archetype scores
remove_ind = which(os_clusters$cluster_annotation %in% c("OC","Myel","TIL","Endo"))
x = as.matrix(arch_osm[,-remove_ind])
rownames(x) = 1:12
cluster_colors = gg_color_hue(nlevels(os_clusters$cluster_annotation)) %>% setNames(levels(os_clusters$cluster_annotation))
clusts = os_clusters$cluster_annotation[-remove_ind] %>% droplevels(c("OC","Myel","TIL","Endo"))
pdf( file.path(outdir, "OS11_Arch_Heatmap.pdf"),
     width=10, height=4)
Heatmap(x, name = "NMF Score", col = col_fun,
        cluster_rows = F, row_names_side = "left",
        row_title = "Archetype", row_title_side = "left",
        column_split = meta_os$tumor_id[-remove_ind], cluster_column_slices = F,
        show_column_names = F, show_column_dend = F,
        column_title_gp = grid::gpar(fontsize = 8), #column_title_rot = 90,
        #top_annotation = HeatmapAnnotation(cell_cluster = clusts,
        #                                   col = list(cell_cluster=cluster_colors))
        )
dev.off()


# figure 4D: OS11 composition plot

# recompute compositions
comp_os2 = data.frame(Archetype = apply(arch_osm[,-remove_ind], 2, which.max),
                      id = meta_os$tumor_id[-remove_ind]) %>%
  mutate(Archetype = factor(Archetype, levels = c(1:12))) %>%
  group_by(id, Archetype) %>% summarize(n = n()) %>%
  reshape2::dcast(Archetype ~ id)
comp_os2[is.na(comp_os2)] = 0
comp_os2[,-1] = lapply(comp_os2[,-1], function(x) x/sum(x))

# save composition barplot with grouped "Other" archetypes: 1,4,7,8,11,12)
pdf( file.path(outdir, "OS11_Arch_Compositions.pdf"),
     width=4, height=4)
# sum "Other" archetypes: 1,4,7,8,11,12
comp_os2 %>% mutate(other = Archetype %in% c(1,4,7,8,11,12)) %>%
  group_by(other) %>% summarise_if(is.numeric,sum) %>% filter(other==T) %>%
  rename(Archetype = other) %>% mutate(Archetype = "Other") %>%
  rbind(comp_os2 %>% filter( !(Archetype %in% c(1,4,7,8,11,12))),.) %>%
  reshape2::melt() %>%
  ggplot() + geom_col(aes(y=variable, x=value, fill = Archetype),
                      position = position_fill(reverse = TRUE)) +
  scale_y_discrete(limits=rev) +
  theme_classic() + ylab("Tumor ID") + xlab("Fraction cells")
dev.off()



