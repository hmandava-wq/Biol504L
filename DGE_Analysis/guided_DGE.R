if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("readr")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("txdbmaker")
BiocManager::install("AnnotationDbi")
BiocManager::install("pheatmap")
BiocManager::install("tximport")


library("tximport")
library("readr")
library("DESeq2")
library("apeglm")
library("txdbmaker")
library("AnnotationDbi")
library("pheatmap")

txdb <- makeTxDbFromGFF(file = "/Users/harikamandava/Desktop/Biol504L/DGE_Analysis/Drosophila.Genome.gtf", format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
head(tx2gene)
write.csv(tx2gene, file = "/Users/harikamandava/Desktop/Biol504L/DGE_Analysis/Drosophila.Genome.gtf", row.names = FALSE)

dir <- "/Users/harikamandava/Desktop/Biol504L/DGE_Analysis/"
samples <- read.csv(file.path(dir, "DGE.Drosophila.Metadata.csv"))

head(samples)
quant_files <- file.path(dir, samples$quant_folder, "quant.sf")
head(quant_files)

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)
samples$group <- factor(paste(samples$Geographic.Location, samples$Temperature, sep = "_"))
samples$group
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)
View(counts(dds))
rld = rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = 'Temperature')
plotPCA(rld, intgroup = c('Temperature', 'Geographic.Location'))
meta <- data.frame(samples$Temperature, samples$Geographic.Location, row.names = colnames(txi$counts))
# converts the dataframe to a matrix
rld_mat <- assay(rld)

# generates the pairwise correlation matrix from the raw matrix
rld_corr <- cor(rld_mat)
pheatmap(rld_corr, annotation = meta)

good_rows <- rowSums(counts(dds)) >= 10
dds <- dds[good_rows,]

dds <- DESeq(dds)
resultsNames(dds)
result_maine_21C_29C <- results(dds, contrast = c("group", "Maine_29C", "Maine_21C"), alpha = 0.05)
summary(result_maine_21C_29C)
res_maine_21_29_apeglm <- lfcShrink(dds, coef = "group_Maine_29C_vs_Maine_21C", type = 'apeglm', res = result_maine_21C_29C)
summary(result_maine_21C_29C)
plotMA(result_maine_21C_29C, ylim=c(-2,2))
plotMA(res_maine_21_29_apeglm, ylim=c(-2,2))
out_csv <- file.path(dir,"Maine.21C_29C.DEG.csv")
write.csv(as.data.frame(res_maine_21_29_apeglm), file = out_csv)