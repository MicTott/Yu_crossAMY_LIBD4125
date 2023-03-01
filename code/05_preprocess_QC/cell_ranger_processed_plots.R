### LS PRELIM snRNA-seq analysis (mouse)
### Build SCE from Cell Ranger-processed count matrices
### qrsh -l bluejay,mf=40G,h_vmem=42G
### Initiated MNT 19Jan2022

# test.edit

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(BiocParallel)
library(scater)
library(scuttle)
library(here)
library(sessioninfo)
library(cowplot)
library(scran)
library(patchwork)

here()
# [1] "/dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125"


## These operations can definitely be built into a function to streamline, but let's
# first read in the count data from the processed (nuclei-already-called) first sample,
# '1M_C_LS':
Sys.time()
# [1] "2022-01-19 13:16:11 EST"
sce.863 <- read10xCounts(samples = "./processed-data/03_cellranger/GSM5836863_aggr/outs/count/filtered_feature_bc_matrix",
                         sample.names = "GSM-863",
                         type = "sparse",
                         col.names = TRUE
                         )


# A summary of this object:
sce.863
# class: SingleCellExperiment 
# dim: 33538 77254 
# metadata(1): Samples
# assays(1): counts
# rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
# ENSG00000268674
# rowData names(3): ID Symbol Type
# colnames(77254): AAACCCAAGAAATGGG-1 AAACCCAAGAGTTGCG-1 ...
# TTTGTTGTCGTAGCTA-3 TTTGTTGTCTGACCCT-3
# colData names(2): Sample Barcode
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):


# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce.863)$Symbol.uniq <- scuttle::uniquifyFeatureNames(ID = rowData(sce.863)$ID, 
                                                              names = rowData(sce.863)$Symbol)
rownames(sce.863) <- rowData(sce.863)$Symbol.uniq


## Read in Cell Ranger's t-SNE and add as an entry to the reducedDims slot, calling it 'TSNE_CR'
TSNE.coords <- read.csv(
  file = "./processed-data/03_cellranger/GSM5836863_aggr/outs//count/analysis/tsne/gene_expression_2_components/projection.csv",
  header = TRUE, row.names = 1
)

head(TSNE.coords)
# TSNE.1    TSNE.2
# AAACCCAAGAAATGGG-1   2.265911  4.977458
# AAACCCAAGAGTTGCG-1 -10.058977 -3.397977
# AAACCCAAGCGATTCT-1 -10.854428 -6.290817
# AAACCCAAGGGATCGT-1   3.573527  3.993352
# AAACCCAAGGGTGAGG-1  26.170249 -9.177906
# AAACCCAAGGTGGGTT-1 -20.510987 -6.811886


# Check it's the same order as the colnames of the SCE
table(rownames(TSNE.coords) == colnames(sce.863)) # all TRUE, good.
reducedDim(sce.863, "TSNE_CR") <- TSNE.coords



## Do the same for the graph-based clustering - in this case as a column vector (class 'factor') in the colData:
graph.clust <- read.csv(
  file = "./processed-data/03_cellranger/GSM5836863_aggr/outs/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
  header = TRUE, row.names = 1
)

# Check it's the same order as the colnames of the SCE
table(rownames(graph.clust) == colnames(sce.863)) # all TRUE, good.
colData(sce.863)$cluster.graph <- factor(graph.clust$Cluster)

# What's the distribution?
table(sce.863$cluster.graph)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 17314  6605  4467  4252  3457  3430  3001  2976  2836  2747  2736  2572  1981 
# 14    15    16    17    18    19    20    21    22    23    24    25    26 
# 1526  1373  1364  1250  1236  1097  1053   990   973   944   889   854   832 
# 27    28    29    30    31    32    33    34 
# 736   676   607   580   544   533   441   382 


## Plot TSNE, coloring by graph-based clusters - this is the same as what's in the web summary :)
# pdf("./plots/02_cellranger_plotting/SRR-863-analysis_TSNE.pdf")
# plotReducedDim(sce.863, dimred="TSNE_CR", colour_by="cluster.graph", text_by="cluster.graph")
# dev.off()

# Can also print some violin plots of your favorite genes - however we first want to log-transform-normalize
#   the counts, due to such variable total counts/cell type(/nuclei)
sce.863 <- logNormCounts(sce.863, assay.type = "counts", log = TRUE, pseudo.count = 1)

pdf("./plots/02_cellranger_plotting/GSM-863-analysis_TSNE_and_Violin_plots.pdf")
# TSNE:
plotReducedDim(sce.863, dimred = "TSNE_CR", colour_by = "cluster.graph", text_by = "cluster.graph") +
  ggtitle("Cell-Ranger-automated analysis: TSNE (sample 'GSM-863')")
# Violin plots of whatever genes you want to print
plotExpression(sce.863,
               x = "cluster.graph", colour_by = "cluster.graph", exprs_values = "logcounts",
               features = c("SNAP25", "SLC17A7", "GAD1", "SST", "CORT", "CRHBP", "NPY"),
               ncol = 2, scales = "free_y"
)
dev.off()


# Using a neater 'plotExpressionCustom' Louise helped me create, for prettier graphics:
#   (italicizes the gene names and plots the median)
source("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/plotExpressionCustom.R")

pdf("./plots/02_cellranger_plotting/custom_broad_markers.pdf", width = 7, height = 10)
plotExpressionCustom(sce.863,
                     anno_name = "cluster.graph", exprs_values = "logcounts",
                     # Broad cell types
                     features = c(
                       "SNAP25", "SYT1", # neuronal
                       "SLC17A7", "SLC7A6", # excit
                       "GAD1", "GAD2", # inhib
                       "SST", "VIP",
                       "CRHBP","LAMP5",
                       "CORT", "PVALB",
                       "NPY", "CCK"
                     ), # mural
                     features_name = "custom-selected",
                     ncol = 2, scales = "free_y"
)
dev.off()


pdf("./plots/02_cellranger_plotting/custom_broad_markers_TSNE.pdf", width = 21, height = 30)
plotlist <- list()
for (i in c("SNAP25", "SYT1", # neuronal
            "SLC17A7", "SLC7A6", # excit
            "GAD1", "GAD2", # inhib
            "SST", "VIP",
            "CRHBP","LAMP5",
            "CORT", "PVALB",
            "NPY", "CCK")) {
  plotlist[[i]] <- plotReducedDim(sce.863, dimred = "TSNE_CR", colour_by = i, by_exprs_values = "logcounts") +
    scale_fill_gradientn(colours = colorRampPalette(c("grey90", "orange3", "firebrick",
                                                              "firebrick", "red", "red"))(10)) + ggtitle(label = i) + 
                                                                theme(plot.title = element_text(size = 20))
}
plot_grid(ncol = 3, plotlist = plotlist)
dev.off()



# ============ Correlation betweeng genes ================

# Computing between specific pairs:
out <- correlatePairs(sce.863, pairings=rbind(c('SST','CORT'), c('SST', 'CRHBP'),c('SST', 'NPY'),c('CORT','CRHBP')), iters=1e5)
head(out)

# Cort scatterplots
pdf("./plots/02_cellranger_plotting/correlate_markers_all_clusters.pdf", width = 7, height = 7)
p1 <- plotExpression(sce.863, features="CORT", x="SST")
p2 <- plotExpression(sce.863, features="CRHBP", x="SST")
p3 <- plotExpression(sce.863, features="NPY", x="SST")
p4 <- plotExpression(sce.863, features="TAC1", x="SST")
p1+p2+p3+p4
dev.off()


# ======== Subset SST ==========
# Computing between specific pairs:
out <- correlatePairs(sce.cluster.32, pairings=rbind(c('SST','CORT'), c('SST', 'CRHBP'),c('SST', 'NPY'),c('CORT','CRHBP')), iters=1e5)
head(out)

# Cort scatterplots
pdf("./plots/02_cellranger_plotting/correlate_markers_cluster32.pdf", width = 7, height = 7)
p1 <- plotExpression(sce.cluster.32, features="CORT", x="SST")
p2 <- plotExpression(sce.cluster.32, features="CRHBP", x="SST")
p3 <- plotExpression(sce.cluster.32, features="NPY", x="SST")
p4 <- plotExpression(sce.cluster.32, features="TAC1", x="SST")
p1+p2+p3+p4
dev.off()
