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

SD_VS_SWW <- results(DESeq_data, 
                     contrast=list(c("treatment_SD_vs_SWW")), 
                     pAdjustMethod="BH")

summary(SD_VS_SWW)
head(SD_VS_SWW)
DESeq2::plotMA(SD_VS_SWW)
DESeq2::plotMA(SD_VS_SWW, main = "DESeq2", ylim=c(-2,2))



SD_VS_SWW<- SD_VS_SWW[order(SD_VS_SWW$padj),]
SD_VS_SWW_dat <- data.frame(SD_VS_SWW@listData$baseMean, SD_VS_SWW@listData$log2FoldChange, SD_VS_SWW@listData$padj )
colnames(SD_VS_SWW_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = SD_VS_SWW$log2FoldChange, x =log2(SD_VS_SWW$baseMean) , 
     ylab = "log2 fold change (SWW/SD)", xlab = "mean expression SWW",
     ylim=c(-2,2), main = "SD differences with SWW", pch = 20, col ="grey")
with(SD_VS_SWW_dat[SD_VS_SWW_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")




SRW_VS_SWW <- results(DESeq_data, 
                    contrast=list(c("treatment_SRW_vs_SWW")), 
                     pAdjustMethod="BH")

summary(SRW_VS_SWW)
DESeq2::plotMA(SRW_VS_SWW)
DESeq2::plotMA(SRW_VS_SWW, main = "DESeq2", ylim=c(-2,2))

SRW_VS_SWW<- SRW_VS_SWW[order(SRW_VS_SWW$padj),]
SRW_VS_SWW_dat <- data.frame(SRW_VS_SWW@listData$baseMean, SRW_VS_SWW@listData$log2FoldChange, SRW_VS_SWW@listData$padj )
colnames(SRW_VS_SWW_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = SRW_VS_SWW$log2FoldChange, x =log2(SRW_VS_SWW$baseMean) , 
     ylab = "log2 fold change (SWW/SRW)", xlab = "mean expression SWW",
     ylim=c(-2,2), main = "SRW differences with SWW", pch = 20, col ="grey")
with(SRW_VS_SWW_dat[SRW_VS_SWW_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")


