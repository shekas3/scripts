library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(gplots)
experimental_design1 <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/experimental_design1.csv", header = TRUE)
modell.counts <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/model1.counts.csv")
rownames(modell.counts) <- modell.counts$X
modell.counts$X <- NULL
experimental_design1$X <- NULL


DESeq_data <- DESeqDataSetFromMatrix(modell.counts, experimental_design1, 
                                     design = formula(~ treatment))
DESeq_data$treatment <- relevel(DESeq_data$treatment, "SWW")
DESeq_data <- DESeq(DESeq_data)
DESeq_data2 <- results(DESeq_data, alpha=0.05)
attr(DESeq_data , "modelMatrixType")
summary(DESeq_data2, alpha= 0.05)
plotDispEsts(DESeq_data)
head(DESeq_data2)
resultsNames(DESeq_data)

res_contrast_simulation_SD_VS_SWW <- results(DESeq_data, contrast=c("treatment", "SD", "SWW"))
summary(res_contrast_simulation_SD_VS_SWW)
head(res_contrast_simulation_SD_VS_SWW)
DESeq2::plotMA(res_contrast_simulation_SD_VS_SWW)


res_contrast_simulation_SRW_VS_SWW <- results(DESeq_data, contrast=c("treatment", "SRW", "SWW"))
summary(res_contrast_simulation_SRW_VS_SWW)
head(res_contrast_simulation_SRW_VS_SWW)
DESeq2::plotMA(res_contrast_simulation_SRW_VS_SWW)

