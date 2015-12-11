library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(gplots)
experimental_design2 <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/experimental_design2.csv", header = TRUE)
model2.counts <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/model2.counts.csv")
rownames(model2.counts) <- model2.counts$X
model2.counts$X <- NULL
experimental_design2$X <- NULL

DESeq_data <- DESeqDataSetFromMatrix(model2.counts, experimental_design2, 
                                     design = formula(~ simulation + treatment + simulation:treatment))
DESeq_data$treatment <- relevel(DESeq_data$treatment, "SWW")
DESeq_data <- DESeq(DESeq_data)
DESeq_data2 <- results(DESeq_data, alpha=0.05)
attr(DESeq_data , "modelMatrixType")
summary(DESeq_data2, alpha= 0.05)
plotDispEsts(DESeq_data)
head(DESeq_data2)
resultsNames(DESeq_data)

res_contrast_simulation_SD_VS_SWW <- results(DESeq_data, 
                                   contrast=list(c("treatment_SD_vs_SWW")), 
                                   pAdjustMethod="BH")
summary(res_contrast_simulation_SD_VS_SWW)
head(res_contrast_simulation_SD_VS_SWW)
DESeq2::plotMA(res_contrast_simulation_SD_VS_SWW)
res_contrast_simulation_SRW_VS_SWW <- results(DESeq_data, 
                                             contrast=list(c("treatment_SRW_vs_SWW")), 
                                             pAdjustMethod="BH")

summary(res_contrast_simulation_SRW_VS_SWW)

DESeq2::plotMA(res_contrast_simulation_SRW_VS_SWW)




