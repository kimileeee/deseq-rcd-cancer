## Downloading RCD signatures from the MSigDB
#. 
#. Access the msigdb documentation here:
#. https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
#. https://igordot.github.io/msigdbr/articles/msigdbr-intro.html

## Pre-requisite: Installing "msigdbr"
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msigdbr")

## Libraries
library(msigdbr)
library(dplyr)

## Downloading data
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
