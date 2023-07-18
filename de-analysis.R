###################
##   Libraries   ##
###################
library(TCGAbiolinks)
library(tidyverse)
library(DESeq2)
library(msigdbr)
library(stringr)


################################
## Downloading TCGA-COAD Data ##
################################
# get a list of projects
#gdcprojects <- getGDCprojects()

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


#############
# Preparing data for DE Analysis #
##########
# prepare count data
tcga_coad_data <- GDCprepare(query_coad, summarizedExperiment = TRUE)
head(tcga_coad_data)

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

##################################################
## Differential Expression Analysis (All Genes) ##
##################################################

# 
dds <- DESeqDataSetFromMatrix(countData = coad_matrix,
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


######################################
## Differential Expression Analysis ##
######################################

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

