# Load packages
library(readr)
library(readxl)
library(dplyr)
library(matrixStats)
library(ggplot2)

# Load Data
X2022_rasagiline <- read_csv("AMED-2022-Rasagiline-Before-After.csv") %>% as.data.frame()
row.names(X2022_rasagiline) <- X2022_rasagiline$Gene 

# Select the 47 PD patients with exosome miRNA data for both before 
# and after rasagiline treatment
X2022_rasagiline <- X2022_rasagiline[,2:95]

# Remove PD10, who did not receive rasailine treatment
X2022_rasagiline <- X2022_rasagiline[,-19:-20]
miRNA <- X2022_rasagiline %>% as.matrix()
row.names(miRNA) <- row.names(X2022_rasagiline)

# miRNA detected in less than half participants of the cohort is regarded as not efficiently detected
miRNA <- miRNA[rowMedians(miRNA) > 0,] 

# For remaining undetected miRNA,
# Their values shall be assigned as half of the minimum value of the detected target miRNA
miRNA[miRNA == 0] <- NA
nas <- which(is.na(miRNA), arr.ind = TRUE)
# 1 is added to avoid the possible production of unwanted 0s
miRNA[nas] <- round((rowMins(miRNA, na.rm = TRUE, useNames = TRUE) [nas[,1]]+1)/2)
miRNA <- miRNA %>% as.data.frame()


# DESeq2 for treatment
library(DESeq2)
library(zoo)
Time <- c(replicate(46,c("Before","After")))
Clinical_data <- read_excel("Clinical-Parameters.xlsx", col_names = FALSE)
Therapy <- Clinical_data[1, 2:95] %>% unlist() %>% na.locf()
Therapy <- Therapy[-19:-20]
Gender <- Clinical_data[2, 2:95] %>% unlist() %>% na.locf()
Gender <- Gender[-19:-20]
HY <- Clinical_data[4,2:95] %>% unlist() %>%  na.locf()
HY <- HY[-19:-20]
HY <- recode(HY, "1" = "I", "2" = "II", "3" = "III")
Age <- Clinical_data[3, 2:95] %>% unlist() %>% as.numeric() %>% na.locf()
for (i in c(1:length(Age))) {
  Age[i] <- if_else( Age[i] < 60, "<60", if_else(Age[i] >70, ">70", "60~70"))
}
Age <- Age[-19:-20]
coldata <- cbind(c(colnames(miRNA)), Time, Gender, HY, Age, Therapy) %>% as.data.frame()
row.names(coldata) <- coldata$V1
colnames(coldata) <- c("Number", "Time", "Gender", "HY", "Age", "Therapy")
dds <- DESeqDataSetFromMatrix(countData = miRNA, colData = coldata, design = ~Time)
dds$Time <- relevel(dds$Time, "Before", "After")
dds <- DESeq(dds)
res <- results(dds)
results <- res %>% as.data.frame()
results <- results %>% 
  mutate(miRNA = row.names(results))
results$miRNA <- gsub(".*_","",results$miRNA)
results$Rank = results$log2FoldChange * abs(log10(results$padj))

# PCA
library(pcaExplorer)
vsdata <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "Time")
plotPCA(vsdata, intgroup = "Gender")
plotPCA(vsdata, intgroup = "HY")
plotPCA(vsdata, intgroup = "Age")
plotPCA(vsdata, intgroup = "Therapy")
pcaExplorer(dds = dds, dst = vsdata)

# DEG
miRNA_Sig_In <- results[results$log2FoldChange > 0 & results$padj <= 0.05, ]
miRNA_Sig_De <- results[results$log2FoldChange < 0 & results$padj <= 0.05, ]
write.csv(miRNA_Sig_In, "Increased-miRNA.csv")
write.csv(miRNA_Sig_De, "Decreased-miRNA.csv")


# VolcanoPlot
library(EnhancedVolcano)
EnhancedVolcano(results, x = "log2FoldChange", y = "padj", lab = results$miRNA,
                pCutoff = 0.05, FCcutoff = 1, legendIconSize = 5, title = NULL,
                subtitle = NULL, gridlines.major = FALSE, gridlines.minor = FALSE, selectLab = FALSE, legendLabels = c("NS", "FC", "pval", "FC&pval"))

write.csv(miRNA, "Countdata_Trimmed.csv")
write.csv(results, "GSEA-Raw-Data.csv")

