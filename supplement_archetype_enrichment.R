# supplemental figure or archetype enrichment analysis
# based on rank correlation of genes with each archetype score (within MTL)


# load in gene correlations (within MTL) for each archetype
# and archetype coefficient weights

library(tidyverse)

setwd('/Users/linclab/Documents/Tannenbaum Lab/SarcomaLandscape')

# specify directory to save enrichment results
resdir = "Results/enrichment"

# load in NMF Archetype scores
datadir = "Results/MTL_Archetype_Results"

# Archetype scores for MTL
mtl_barcodes = read_csv( file.path(datadir, "MTL_Barcodes.csv"), col_names = F) %>% pull(1)
arch_mtl = read_csv( file.path(datadir, "MTL_Arch.csv"), col_names = F) %>% as.matrix
rownames(arch_mtl) = paste0("Arch",1:12)
colnames(arch_mtl) = mtl_barcodes
nArch = dim(arch_mtl)[1]

# cell metadata for each dataset
# match MTL metadata to barcodes used in archetype analysis
meta_mtl = read_delim("Data/MTL_OAC_SCT_qnorm/cell_metadata_clustered.txt")
meta_mtl_matched = meta_mtl[match(mtl_barcodes, meta_mtl$barcode),] %>%
  mutate(cluster_annotation = factor(cluster_annotation,
                                     levels = c("UD", "MSC-H", "MSC-L",
                                                "O1", "O2",
                                                "A1", "A2", "A3", "A4",
                                                "CP", "C1", "C2", "C3")))

# load gene expression data from MTL (combined with chondro, SCT normalized)
library(Seurat)
MTL_obj = readRDS( "Data/MTL_OAC_SCT_qnorm/MTL_OAC_harmony.rds" ) %>%
  subset(cells = mtl_barcodes) # pull only relevant cells (n=22589 cells)


# load archetype coefficient weights
arch_scores = read_csv( file.path(datadir, "Archetype_Gene_Scores.csv"))
arch_scores_mat = as.matrix(arch_scores[,-1])
rownames(arch_scores_mat) = arch_scores$Gene


# heatmap of NMF loadings
col_fun = colorRamp2(c(0, 0.01, 0.02), c("#009FFF", "#FFFF00", "#FF5100"))
vec = 1:12; wtavg = rowSums(sweep(arch_scores_mat, MARGIN=2, vec, `*`)) / rowSums(arch_scores_mat)
pdf(file.path(resdir,'archetype_loadings_heatmap.pdf'), height=20, width=5)
Heatmap(arch_scores_mat, name="NMF loading", col=col_fun,
        cluster_columns = F, column_names_side = "top",
        column_title = "Archetype", column_names_rot = 0,
        cluster_rows = T, row_dend_reorder = wtavg, clustering_method_rows = "ward.D2",
        row_names_gp = gpar(fontsize=3))
dev.off()



##### CORRELATION ANALYSIS #####

# recompute correlation of all genes with each archetype
expr = as.matrix(MTL_obj@assays$SCT@data) # using SCTransform normalized expression
expr = expr[rowSums(expr) > 0,] # remove genes with 0 expression
dim(expr) # 27639 x 22589

# compute Gene x Archetype correlation
# here considering both Pearson and Spearman
t0 = proc.time()
r_p = cor( t(expr), t(arch_mtl), method = "pearson" )
rho = cor( t(expr), t(arch_mtl), method = "spearman" )
proc.time() - t0 # takes about 2-3 min


thresh = apply(rho, 2, quantile, probs = c(0.95),  na.rm = TRUE)



# heatmap of Spearman correlation
library(ComplexHeatmap)
library(circlize)
pdf('test.pdf', width=20, height=5)
Heatmap(t(rho), name="Spearman rho",
        cluster_rows = F, row_names_side = "left",
        row_title = "Archetype", row_title_side = "left",
        show_column_names = F,
        cluster_columns = T)#, column_names_gp = gpar(fontsize=3))
dev.off()


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

rho_subset = rho[selected_genes,]
r_subset = r_p[selected_genes,]

# heatmaps for both pearson and spearman correlation of archetype x gene expression
pdf( file.path(resdir, "archcorr_spearman_subset.pdf"), width=6, height=10)
Heatmap(rho_subset, name="Spearman corr.",
        row_split = gene_class, cluster_row_slices = F, # cluster_rows = F,
        row_title_rot = 30, row_title_side = "right", show_row_dend = F, cluster_rows=F,
        cluster_columns = F,  column_names_side = "top")
dev.off()
pdf( file.path(resdir, "archcorr_pearson_subset.pdf"), width=6, height=10)
Heatmap(r_subset, name="Pearson corr.",
        row_split = gene_class, cluster_row_slices = F, # cluster_rows = F,
        row_title_rot = 30, row_title_side = "right", show_row_dend = F, cluster_rows=F,
        cluster_columns = F, column_names_side = "top")
dev.off()


# another time for archetype loadings
mat2 = arch_scores_mat[selected_genes[which(selected_genes %in% rownames(arch_scores_mat))],]
gene_class2 = gene_class[which(selected_genes %in% rownames(arch_scores_mat))]
pdf( file.path(resdir, "archetype_loadings_subset.pdf"), width=6, height=10)
Heatmap(mat2, name="NMF loading", col=col_fun,
        cluster_columns = F, column_names_side = "top",
        column_title = "Archetype", column_names_rot = 0, column_names_centered = T,
        row_names_side = "right", row_title_side = "right",
        row_split = gene_class2, row_title_rot = 45,
        cluster_rows = F, cluster_row_slices = F)
dev.off()


# another time for archetype loadings
mat3 = arch_scores_mat[selected_genes[which(selected_genes %in% rownames(arch_scores_mat))],]
gene_class3 = gene_class[which(selected_genes %in% rownames(arch_scores_mat))]
pdf( file.path(resdir, "archetype_loadings_subset.pdf"), width=6, height=10)
Heatmap(mat2, name="NMF loading", col=col_fun,
        cluster_columns = F, column_names_side = "top",
        column_title = "Archetype", column_names_rot = 0, column_names_centered = T,
        row_names_side = "right", row_title_side = "right",
        row_split = gene_class2, row_title_rot = 45,
        cluster_rows = F, cluster_row_slices = F)
dev.off()


# we convert spearman rho into a p-value using t statistic
# reference: https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
#  Kendall, M. G.; Stuart, A. (1973). "Sections 31.19, 31.21". The Advanced Theory of Statistics, Volume 2: Inference and Relationship. Griffin. ISBN 978-0-85264-215-3.
n = dim(expr)[2] # sample size = number of cells = 22589

t_p = r_p * sqrt( (n-2)/(1-r_p^2) ) # t-statistic for Pearson correlation
t_s = rho * sqrt( (n-2)/(1-rho^2) ) # t-statistic for Spearman correlation
p_p = pt( t_p, n-2, lower.tail = F) # here, we consider a right sided test because we want positive correlations
p_s = pt( t_s, n-2, lower.tail = F) # here, we consider a right sided test because we want positive correlations

padj_p = p_p %>% as.data.frame %>%
  mutate(across(contains("Arch"), ~ p.adjust(.x, method="BH", n=length(.x))))
padj_s = p_s %>% as.data.frame %>%
  mutate(across(contains("Arch"), ~ p.adjust(.x, method="BH", n=length(.x))))

colSums(p_p<0.05)
colSums(padj_p<0.05)
colSums(p_p<0.01)
colSums(padj_p<0.01)
colSums(p_p<0.001)
colSums(padj_p<0.001)

colSums(p_s<0.05)
colSums(padj_s<0.05)
colSums(p_s<0.01)
colSums(padj_s<0.01)
colSums(p_s<0.001)
colSums(padj_s<0.001)


# looking at a few correlated example pairs
r_p[8:12,14] # not much correlation
r_p[8:12,15] # good correlation
cor.test( as.numeric(arch_mtl[8,]), as.numeric(expr[15,]) , alternative = "g")
cor.test( as.numeric(arch_mtl[9,]), as.numeric(expr[15,]) )
cor.test( as.numeric(arch_mtl[10,]), as.numeric(expr[15,]), alternative = "g" )
cor.test( as.numeric(arch_mtl[11,]), as.numeric(expr[15,]) )
cor.test( as.numeric(arch_mtl[12,]), as.numeric(expr[15,]) )
(1-p_p[8,14])*2
r_p[8,14]
p_p[8,14]
r_p[8,15]
p_p[10,15]


data.frame(r = rho[,1], pval = p_s[,1],
           padj = padj_s[,1]) %>%
  ggplot() +
  geom_point(aes(x=r,y=pval))
  #geom_histogram(aes(x=r))
data.frame(r = r_p[,1], pval = p_p[,1],
           padj = padj_p[,1]) %>%
  ggplot() +
  geom_point(aes(x=r,y= log10(padj) ))


# 1. GSEA
# determine ranked list for each archetype
# use R package fgsea, msigdbr
library(fgsea)
library(msigdbr)

# load MSigDB hallmark pathway collection
msigdbr_df = msigdbr(species = "human", category = "H") # hallmark
msigdbr_df = msigdbr(species = "human", category = "C2") # curated
msigdbr_df = msigdbr(species = "human", category = "C5") # GO
msigdbr_df = msigdbr(species = "human", category = "C8") # cell type markers
msigdbr_df$gs_subcat %>% unique
# look into if can filter down by subcategory or by differentiation pathways, etc
msigdbr_df = msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS") # curated sub categories: CGP CP:BIOCARTA CP:KEGG CP CP:PID CP:REACTOME CP:WIKIPATHWAYS"

pathways = split(x = msigdbr_df$human_gene_symbol, f = msigdbr_df$gs_name)


# use correlation as ranked statistic
res_list_gsea = list()
plot_list_gsea = list()
for (iArch in 1:nArch) {
  rankData = rho[,iArch] %>% sort(decreasing = T)
  
  # run fgsea
  res = fgsea(pathways,
               rankData,
               minSize = 15,
               maxSize = 500,
               nperm = 1000)
  #res %>% arrange(pval) %>% pull(pathway) %>% head
  #res$pathway[res$padj < 0.05] %>% head
  print(sprintf("Arch %i: %i sig. pathways (FDR < 0.05)", iArch, sum(res$padj<0.05,na.rm=T)))
  #print( res$pathway[res$padj < 0.05] %>% paste(collapse = "; ") )
  
  res_list_gsea[[iArch]] = res
  topUp <- res %>%
    filter(ES > 0) %>%
    top_n(10, wt=-padj) %>%
    top_n(10, wt=ES) %>%
    arrange(-ES)
  plot_list_gsea[[iArch]] = plotGseaTable(pathways[topUp$pathway],
                                          rankData,
                                          res,
                                          gseaParam = 0.5,
                                          render=F)
}
p_enrich = ggpubr::ggarrange(plotlist=plot_list_gsea, nrow=3, ncol=4)
p_enrich


# save each archetype enrichment to a excel worksheet page
library(xlsx)
filename = file.path(resdir,"archetype_GSEA.xlsx")
write.xlsx(df1, file=filename, sheetName="Arch1_GSEA_H", row.names=FALSE)
for (i in 2:narch) {
  write.xlsx(df2, file=filename,
             sheetName=paste0("Arch",i,"_GSEA_H"), row.names=FALSE,
             append=TRUE)
}


# 2. GO enrichment (using clusterProfiler package)
# based on list of significantly correlated genes (by Spearman test) to each archetype
library(clusterProfiler)
library(enrichplot)

# first convert HGNC gene symbol to Entrez gene ID
library(org.Hs.eg.db)
#genes_HGNC = rownames(expr) # all gene expression
genes_HGNC = rownames(arch_scores_mat) # archetype loadings
genes_HGNC = rownames(rho) # spearman correlation
genes_entrez = mapIds(org.Hs.eg.db,
                      keys = genes_HGNC,
                      keytype = "SYMBOL",
                      column = "ENTREZID",
                      multiVals = "first")
#arch_scores_sub = arch_scores_mat[!is.na(genes_entrez),]
#genes_entrez_sub = genes_entrez[!is.na(genes_entrez)]

# compute GO enrichment of each archetype
# using a cutoff of adjusted p-value
# or top percentile of correlated genes
# or based on archetype coefficients

ego_list = list() # store each archetype enrichment results in list
plot_list = list() # create enrichment dotplots of each archetype
arch_gene_list = list() # store significant genes
for (iArch in 1:nArch) {
  #arch_genes = genes_entrez[padj_s[,iArch] < 0.001] %>% na.omit # entrezid of correlated genes
  #arch_genes = data.frame(r = rho[,iArch], entrez = genes_entrez) %>% na.omit %>%
  #  arrange(-r) %>% top_frac(0.05, r) %>% pull(entrez) # take top 1% correlated genes
  #arch_genes = data.frame(c = arch_scores_mat[,iArch], entrez = genes_entrez) %>% na.omit %>%
  #  arrange(-c) %>% top_frac(0.05, c) %>% pull(entrez) # take 5% top archetype scores
  
  #thresh = quantile(arch_scores_sub[,iArch],0.8) # take top 25 percentile of archetype loadings
  #arch_genes = genes_entrez_sub[which(arch_scores_sub[,iArch] > thresh)]
  
  #thresh = quantile(r_p[,iArch],0.95) # take top 25 percentile of archetype loadings
  thresh = 0.2
  arch_genes = genes_entrez[which(rho[,iArch] > thresh)]
  
  print(paste("Archetype",iArch,": ngenes =",length(arch_genes)))
  ego = enrichGO( gene = arch_genes %>% na.omit,
                  universe = genes_entrez %>% na.omit,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = T)
  ego_list[[iArch]] = ego
  print( paste(" Arch",iArch,": n pathways =", dim(ego)[1]))
}
# create enrichment dotplots of each cluster
for (iArch in 1:nArch) {
  thresh = 0.2
  arch_genes = genes_entrez[which(rho[,iArch] > thresh)]
  arch_gene_list[[iArch]] = arch_genes
  
  if (dim(ego_list[[iArch]])[1]==0) {plot_list[[iArch]] = patchwork::plot_spacer(); next} # if no pathways, plot a spacer
  plot_list[[iArch]] = dotplot(ego_list[[iArch]], font.size = 8, showCategory=10 ) +
    ggtitle(sprintf( "Archetype %i (n=%i)", iArch, length(arch_genes) )) + # length(genes_entrez[padj_s[,iArch] < 0.001] %>% na.omit)
    scale_color_gradient(low="red", high="blue", limits = c(NA, 0.01),
                         trans =  "log10",
                         labels = function(x) sprintf("%.1e",x)) +
    theme(axis.text.y = element_text(size=6),
          plot.title = element_text(size = 10))
}
# combine into one plot
p_enrich = ggpubr::ggarrange(plotlist=plot_list, nrow=3, ncol=4,
                             common.legend=T, legend="right")
#pdf( file.path(resdir, "MTL_Archetype_EnrichGO.pdf"), width=20, height=10)
p_enrich
#dev.off()


# save each archetype enrichment to a excel worksheet page
library(xlsx)
filename = file.path(resdir,"MTLarchetype_enrichGO.xlsx")
for (iArch in 1:nArch) {
  if (iArch == 1) {append_page = F} else {append_page = T}
  ego_list[[iArch]] %>% as.data.frame %>%
    write.xlsx(file=filename,
               sheetName=paste0("Arch",iArch,"_enrichGO"), row.names=FALSE,
               append=append_page)
}


