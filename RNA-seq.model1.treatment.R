library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(gplots)

experimental_design1 <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/experimental_design1.csv", header = TRUE)
modell.counts <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/model1.counts.csv")
rownames(modell.counts) <- modell.counts$X
modell.counts$X <- NULL
experimental_design1$X <- NULL

#Design DEseq formula based on treatments, ignoring lane effect
DESeq_data <- DESeqDataSetFromMatrix(modell.counts, experimental_design1, 
                                     design = formula(~ treatment))
DESeq_data$treatment <- relevel(DESeq_data$treatment, "SWW") # set SWW treatmmnet as reference treatment ( control)
DESeq_data <- DESeq(DESeq_data) # Clculate differential gene expression 
DESeq_data2 <- results(DESeq_data, alpha=0.05) # incorporate p-values less that 0.05
attr(DESeq_data , "modelMatrixType")
summary(DESeq_data2, alpha= 0.05) # summary of results
plotDispEsts(DESeq_data) # plot the results
head(DESeq_data2)
resultsNames(DESeq_data)

SD_VS_SWW <- results(DESeq_data, contrast=c("treatment", "SD", "SWW")) # DE specifically between SD and SWW treatment
summary(SD_VS_SWW)
head(SD_VS_SWW)
DESeq2::plotMA(SD_VS_SWW) # plot the contrast results

# make a custome data.frame to create custom MAplot which is more useful than plotMA

SD_VS_SWW<- SD_VS_SWW[order(SD_VS_SWW$padj),]
SD_VS_SWW_dat <- data.frame(SD_VS_SWW@listData$baseMean, SD_VS_SWW@listData$log2FoldChange, SD_VS_SWW@listData$padj )
colnames(SD_VS_SWW_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = SD_VS_SWW$log2FoldChange, x =log2(SD_VS_SWW$baseMean) , 
     ylab = "log2 fold change (SWW/SD)", xlab = "mean expression SWW",
     main = "SD differences with SWW", pch = 20, col ="grey")
with(SD_VS_SWW_dat[SD_VS_SWW_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")


SRW_VS_SWW <- results(DESeq_data, contrast=c("treatment", "SRW", "SWW")) # DE specifically between SRW and SWW treatment
summary(SRW_VS_SWW)
head(SRW_VS_SWW)
DESeq2::plotMA(SRW_VS_SWW) # plot the contrast results

# make a custome data.frame to create custom MAplot which is more useful than plotMA
SRW_VS_SWW<- SRW_VS_SWW[order(SRW_VS_SWW$padj),]
SRW_VS_SWW_dat <- data.frame(SRW_VS_SWW@listData$baseMean, SRW_VS_SWW@listData$log2FoldChange, SRW_VS_SWW@listData$padj )
colnames(SRW_VS_SWW_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = SRW_VS_SWW$log2FoldChange, x =log2(SRW_VS_SWW$baseMean) , 
     ylab = "log2 fold change (SWW/SRW)", xlab = "mean expression SWW",
     main = "SRW differences with SWW", pch = 20, col ="grey")
with(SRW_VS_SWW_dat[SRW_VS_SWW_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")





