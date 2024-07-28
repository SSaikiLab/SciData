# Load library
library(readr)
library(readxl)
library(dplyr)
library(matrixStats)
library(ggplot2)

# Load Data
miRNA_2021 <- read_excel("2021.xlsx") %>% as.data.frame()
colnames(miRNA_2021) <- c("miRNA", c(1:76))
row.names(miRNA_2021) <- miRNA_2021$miRNA
miRNA_2021 <- miRNA_2021[-1,-1]


# Remove miRNA detected in less than half participants of the cohort
min(miRNA_2021[miRNA_2021 > 0])
miRNA <- apply(miRNA_2021, 2, as.numeric)
miRNA <- apply(miRNA, 2, as.integer)
row.names(miRNA) <- row.names(miRNA_2021)

# Remove miRNA detected in less than half of the participants  cohort
miRNA <- miRNA[rowMedians(miRNA) > 0,] 

# For remaining undetected miRNA,
# Their values shall be assigned as half of the minimum value of the detected target miRNA
miRNA[miRNA == 0] <- NA
nas <- which(is.na(miRNA), arr.ind = TRUE)
miRNA[nas] <- round((rowMins(miRNA, na.rm = TRUE, useNames = TRUE) [nas[,1]]+1)/2)
miRNA <- miRNA %>% as.data.frame()


# DESeq2 
library(DESeq2)
coldata <- cbind(colnames(miRNA), c(replicate(46, "PD"), replicate(30, "HC")))
colnames(coldata) <- c("Patient", "Disease")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = miRNA, colData = coldata, design = ~Disease)
dds$Disease <- relevel(dds$Disease, "HC", "PD")
dds <- DESeq2::DESeq(dds)
res <- results(dds)
results <- res %>% as.data.frame()
results <- results %>% 
  mutate(miRNA = row.names(results))
results$miRNA <- gsub(".*_","",results$miRNA)

# PCA
library(pcaExplorer)
vsdata <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "Disease")
pcaExplorer(dds = dds, dst = vsdata)


# Volcanoplot
library(EnhancedVolcano)
EnhancedVolcano(results, x = "log2FoldChange", y = "padj", lab = results$miRNA, selectLab = FALSE)

# DEG
miRNA_Sig_In <- results[results$log2FoldChange > 0 & results$padj <= 0.05, ]
miRNA_Sig_De <- results[results$log2FoldChange < 0 & results$padj <= 0.05, ]
write.csv(miRNA_Sig_In, "Increased-miRNA.csv")
write.csv(miRNA_Sig_De, "Decreased-miRNA.csv")

results$Rank <- results$log2FoldChange * (-log10(results$padj))
write.csv(results, "results.csv")
