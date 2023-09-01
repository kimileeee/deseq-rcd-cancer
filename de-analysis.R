#######################
##   Intstallation   ##
#######################
BiocManager::install("TCGAbiolinks")
BiocManager::install("msigdbr")
BiocManager::install("RNAseqQC")
BiocManager::install("DESeq2")
BiocManager::install("ensembldb")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("purrr")
install.packages("tidyr")
install.packages("tibble")
install.packages("magrittr")
install.packages("stringr")
install.packages("tidyverse")
install.packages("vsn")

###################
##   Libraries   ##
###################
library(TCGAbiolinks)
library(msigdbr)
library(RNAseqQC)
library(DESeq2)
library(ensembldb)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(tibble)
library(magrittr)
library(stringr)
library(tidyverse)
library(vsn)

################################
## Downloading TCGA-COAD Data ##
################################
# get a list of projects
# gdcprojects <- getGDCprojects()

# Querying All Colorectal Cancer Data (Tumor)
query_tumor <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open",
  sample.type = "Primary Tumor"
)
tumor <- getResults(query_tumor)

# Querying All Colorectal Cancer Data (Normal)
query_normal <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open",
  sample.type = "Solid Tissue Normal"
)
normal <- getResults(query_normal)

# Filtering samples with both tumor and normal tissues
submitter_ids <- inner_join(tumor, normal, by = "cases.submitter_id") %>%
              select(cases.submitter_id)
tumor <- tumor %>%
  filter(cases.submitter_id %in% submitter_ids$cases.submitter_id)
normal <- normal %>%
  filter(cases.submitter_id %in% submitter_ids$cases.submitter_id)

# Merging to a single data.frame
samples <- rbind(tumor, normal)
unique(samples$sample_type)
samples$sample.submitter_id

# Overwrite barcode of query with the filtered sample ids
query_coad <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open",
  sample.type = c("Solid Tissue Normal", "Primary Tumor"),
  barcode = as.list(samples$sample.submitter_id)
)
getResults(query_coad)

# Downloading TCGA data
GDCdownload(query_coad)


#####################################
## Downloading RCD signatures data ##
#####################################
# Necroptosis
necroptosis_geneset <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")  %>%
  filter(gs_name == "GOBP_NECROPTOTIC_SIGNALING_PATHWAY")
head(necroptosis_geneset)

# Ferroptosis
ferroptosis_geneset <- msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS")  %>%
  filter(gs_name == "WP_FERROPTOSIS")
head(ferroptosis_geneset)

# Ferroptosis
pyroptosis_geneset <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")  %>%
  filter(gs_name == "REACTOME_PYROPTOSIS")
head(pyroptosis_geneset)


#####################################
## Data wranggling for DE Analysis ##
#####################################
# prepare count data
tcga_coad_data <- GDCprepare(query_coad, summarizedExperiment = TRUE)
head(tcga_coad_data)
#tcga_coad_data <- tcga_coad_data %>%
#                    filter()

# unstanded, stranded_first, stranded_second, tpm_unstrand, fpkm_unstrand, fpkm_uq_unstrand
coad_matrix <- assay(tcga_coad_data, 'unstranded')
coad_matrix

# make dataframe of sample names, cell Lines, and tissue type(normal or tumor)
rownames(samples) <- samples$cases
View(samples)
samples <- samples %>%
  select(case="cases.submitter_id", type="sample_type")
samples$type <- str_replace(samples$type, "Solid Tissue Normal", "normal")
samples$type <- str_replace(samples$type, "Primary Tumor", "tumor")

coad_matrix <- coad_matrix[, rownames(samples)]
# Check if all samples in the counts dataframe are in the samples dataframe
all(colnames(coad_matrix) %in% rownames((samples)))
all(colnames(coad_matrix) == rownames(samples))

a <- str_split_fixed(rownames(coad_matrix), "[.]", 2)
counts_matrix <- coad_matrix
rownames(counts_matrix) <- a[, 1]
counts_matrix

#####################
## Quality Control ##
#####################
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = samples,
                              design = ~ type)

# QC plots on raw count matrix
plot_total_counts(dds)
plot_library_complexity(dds)
plot_gene_detection(dds)

# Gene biotypes
plot_biotypes(dds)

# Gene filtering
# filter to keep only rows that have 10 reads total. before: 45194 x 87 after: 53996 x 519
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# dds2 <- filter_genes(dds, min_count = 10, min_rep = 41)

# Variance stabilization
vsd <-vst(dds)
sdplot <-meanSdPlot(assay(vsd))
sdplot$gg <- sdplot$gg +
  ggtitle(label="Mean SD of transformed read counts") +
  ylab("standard deviation")
print(sdplot$gg)

# Chromosomal expression
map(c("1", "5", "14"), ~plot_chromosome(vsd, .x))

# Replicate variability
# define new grouping variable
colData(vsd)$trt_mut <- paste0(colData(vsd)$type, "_", colData(vsd)$case)

ma_plots <- plot_sample_MAs(vsd, group = "trt_mut")
cowplot::plot_grid(plotlist = ma_plots[17:24], ncol = 2)

# Clustering
# set seed to control random annotation colors
set.seed(1)
plot_sample_clustering(vsd, anno_vars = c("type"), distance = "euclidean")


# Principal component analysis (PCA) plot (from Tutorial)
plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "type")

# PCA plot (from Ms Jenny)
pcaData<-plotPCA(vsd, intgroup=c("type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


#to further check outliers
boxplot(log10(assays(dds$type=="normal")[["cooks"]]), range=0, las=2)


##################################################
## Differential Expression Analysis (All Genes) ##
##################################################
library(EnsDb.Hsapiens.v79)
# Plot a gene
dds <- estimateSizeFactors(dds)
# dds1 <- dds
# geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(dds), keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
# geneIDs
# rownames(dds1)
# rownames(dds1) <- geneIDs$SYMBOL[rownames(dds1) == geneIDs$GENEID]
# "TLR3" %in% rownames(dds1)
# dds
# plot_gene("RIPK3", dds1, x_var = "case", color_by = "type")

# Differential expression testing
dds$type <- factor(dds$type, levels = c("normal","tumor")) 
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
res
write.csv(res, file = "DE_Results.csv")

plotDispEsts(dds)

# Plot a testing result
de_res <- lfcShrink(dds, coef="type_tumor_vs_normal", lfcThreshold = log2(1.5), type = "normal", parallel = TRUE)


# https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
res <- results(dds, contrast=c('type', 'tumor', 'normal'))
res <- res[order(res$padj),]
library(knitr)
kable(res[1:5,-(3:4)])

res <- results(dds, contrast=c("type","normal","tumor"))
ix = which.min(res$padj)
res <- res[order(res$padj),]
kable(res[1:5,-(3:4)])

barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )



#
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

resultsNames(dds)

# Visualization MA plot
plotMA(res)


############################################
# Filtering Necroptosis-related regulators #
############################################
a <- str_split_fixed(rownames(coad_matrix), "[.]", 2)
counts_matrix <- coad_matrix
rownames(counts_matrix) <- a[, 1]
counts_matrix

coad_necroptosis <- counts_matrix[rownames(counts_matrix) %in% necroptosis_geneset$ensembl_gene, ]
coad_necroptosis <- coad_necroptosis[, rownames(samples)]

# Check if all samples in the counts dataframe are in the samples dataframe
all(colnames(coad_necroptosis) %in% rownames((samples)))
all(colnames(coad_necroptosis) == rownames(samples))


###############################
## DE Analysis (Necroptosis) ##
###############################

# 
dds <- DESeqDataSetFromMatrix(countData = coad_necroptosis,
                              colData = samples,
                              design = ~ type)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# 
dds$type <- relevel(dds$type, ref = "normal")

#
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

resultsNames(dds)
# Code below allows for comparing across more than 1 levels
# results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))
res_multiomics <- results(dds, lfcThreshold= 2, alpha=0.05)
summary(res_multiomics)

# Visualization MA plot
plotMA(res)

