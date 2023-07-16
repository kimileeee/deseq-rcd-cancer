library(airway)

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

counts_data_sample <- read.csv("counts_data.csv")
head(counts_data_sample)

colData <- read.csv("sample_info.csv")

all(colnames(counts_data_sample) %in% rownames((colData)))

all(colnames(counts_data_sample) == rownames(colData))

dds <- DESeqDataSetFromMatrix(countData = counts_data_sample,
                        colData = colData,
                        design = ~ dexamethasone)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")


dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

resultsNames(dds)
# Code below allows for comparing across more than 1 levels
results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# Visualization MA plot
plotMA(res)
