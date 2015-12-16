library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(gplots)
#Importing the counts from model1 and experimental model 1
experimental_design1 <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/experimental_design1.csv", header = TRUE)
modell.counts <- read.csv("c:/Users/guest6/Desktop/FinalProject/results/model1.counts.csv")
rownames(modell.counts) <- modell.counts$X
modell.counts$X <- NULL
experimental_design1$X <- NULL
#Input the data into DESeq package: it starts with formula
#first we check for lane effect
test_lane_effects <- DESeqDataSetFromMatrix(modell.counts, experimental_design1, 
                                            design = formula(~ lane))

test_lane_effects2 <- DESeq(test_lane_effects) # calculate differential gene expression
test_lane_effects2_results <- results(test_lane_effects2) # results
#up regulated and down regulated log2 fold changes of genes from lane to lane
summary(test_lane_effects2_results) # summary of results 
plotDispEsts(test_lane_effects2) 
#Principal Components analysis and hierarchical clustering to visualize patterns (and to identify potential confounds)
for_pca <- rlog(test_lane_effects2, blind=TRUE) 
plotPCA(for_pca, intgroup=c("lane")) # slightly lane effect

plotPCA(for_pca, intgroup=c("replicate", "barcode"))

rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion
distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix

rownames(mat) <- colnames(mat) <- with(colData(test_lane_effects2), paste(replicate, barcode, lane, sep=" : "))


hc <- hclust(distsRL)  # performs hierarchical clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)  # picking our colours

#Now we generate the plot
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(17, 17))
