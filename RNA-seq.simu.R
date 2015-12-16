library(splitstackshape)
library(dplyr)
setwd("~/bio720/counts_htseq/")
path = "~/bio720/counts_htseq/"
# import the data from HTseq
# create a forlooop to tidy up the data and perform simulation
file.names <- dir(path, pattern = ".htsq")
file.counts <- list ()
total.counts <- list ()
unique.mapped <- list ()
resample1 <- list()
resample2 <- list()
resample3 <- list()
resample4 <- list()
resample5 <- list()
simuSamples1 <- list()
simuSamples2 <- list()
simuSamples3 <- list()
simuSamples4 <- list()
simuSamples5 <- list()
for(i in 1:length(file.names)){
  file.counts[[ file.names[i] ]] <- read.table(file.names[i], header = FALSE, sep = "\t", col.names = c("genes", "count")) # read the data
  total.counts[[ file.names[i] ]] <- expandRows(file.counts[[i]], "count") # this function from splitstackshape replicate a gene name based on it's count
  unique.mapped[[ file.names[i] ]] <- nrow(total.counts[[i]]) # calculate total number of counts in each sample for re-sampling
  resample1[[ file.names[i] ]] <- total.counts[[i]][sample(unique.mapped[[i]], replace = TRUE), ] #re-sampling from total counts without replacement (Bootstrap)
  resample2[[ file.names[i] ]] <- total.counts[[i]][sample(unique.mapped[[i]], replace = TRUE), ] #Five round of re-sampling
  resample3[[ file.names[i] ]] <- total.counts[[i]][sample(unique.mapped[[i]], replace = TRUE), ]
  resample4[[ file.names[i] ]] <- total.counts[[i]][sample(unique.mapped[[i]], replace = TRUE), ]
  resample5[[ file.names[i] ]] <- total.counts[[i]][sample(unique.mapped[[i]], replace = TRUE), ]
  simuSamples1[[ file.names[i] ]] <- data.frame(table(resample1[[i]])) # reshape the resample data based on unique gene names
  simuSamples2[[ file.names[i] ]] <- data.frame(table(resample2[[i]]))
  simuSamples3[[ file.names[i] ]] <- data.frame(table(resample3[[i]]))
  simuSamples4[[ file.names[i] ]] <- data.frame(table(resample4[[i]]))
  simuSamples5[[ file.names[i] ]] <- data.frame(table(resample4[[i]]))
  colnames(simuSamples1[[file.names[i]]]) <- c("genes", "count")
  colnames(simuSamples2[[file.names[i]]]) <- c("genes", "count")
  colnames(simuSamples3[[file.names[i]]]) <- c("genes", "count")
  colnames(simuSamples4[[file.names[i]]]) <- c("genes", "count")
  colnames(simuSamples5[[file.names[i]]]) <- c("genes", "count")
}
#unlist and combination of count data
tot_count_sample <- matrix(unlist(lapply(file.counts, function(x) x$count)), nrow=26351, ncol=6)
colnames(tot_count_sample) <- c("SD_1_R", "SD_2_R", "SRW_1_R", "SRW_2_R", "SWW_1_R", "SWW_2_R")
#total gene names from htseq count
total_genes <- matrix(unlist(lapply(file.counts, function(x) x$gene)), nrow=26351, ncol=1)
colnames(total_genes) <- "genes"
#combination of count data from First round of simulation
tot_count_simu1 <- matrix(unlist(lapply(simuSamples1, function(x) x$count)), nrow=26351, ncol=6)
colnames(tot_count_simu1) <- c("SD_1_S1", "SD_2_S1", "SRW_1_S1", "SRW_2_S1", "SWW_1_S1", "SWW_2_S1")
#combination of count data from Second round of simulation
tot_count_simu2 <- matrix(unlist(lapply(simuSamples2, function(x) x$count)), nrow=26351, ncol=6)
colnames(tot_count_simu2) <- c("SD_1_S2", "SD_2_S2", "SRW_1_S2", "SRW_2_S2", "SWW_1_S2", "SWW_2_S2")
#combination of count data from Third round of simulation
tot_count_simu3 <- matrix(unlist(lapply(simuSamples3, function(x) x$count)), nrow=26351, ncol=6)
colnames(tot_count_simu3) <- c("SD_1_S3", "SD_2_S3", "SRW_1_S3", "SRW_2_S3", "SWW_1_S3", "SWW_2_S3")
#combination of count data from Third round of simulation
tot_count_simu4 <- matrix(unlist(lapply(simuSamples4, function(x) x$count)), nrow=26351, ncol=6)
colnames(tot_count_simu4) <- c("SD_1_S4", "SD_2_S4", "SRW_1_S4", "SRW_2_S4", "SWW_1_S4", "SWW_2_S4")
#combination of count data from Third round of simulation
tot_count_simu5 <- matrix(unlist(lapply(simuSamples5, function(x) x$count)), nrow=26351, ncol=6)
colnames(tot_count_simu5) <- c("SD_1_S5", "SD_2_S5", "SRW_1_S5", "SRW_2_S5", "SWW_1_S5", "SWW_2_S5")
######################################################################################################################
# First MODEL : real samples : Test lane effect
######################################################################################################################
model1.counts <-  data.frame(tot_count_sample)
#Setting up the experimental design object
parse_names <- strsplit(file.names, split="_")
#a matrix of the names of treatment groups
parse_names <- matrix(unlist(parse_names), nrow=6, ncol=5, byrow=T)
col_names_counts <- paste(parse_names[,1], "_", parse_names[,2], "_", parse_names[,3], "_", parse_names[,4], sep="")
colnames(model1.counts) = col_names_counts # sample names as column names
# add gene names as row names
rownames(model1.counts) <- total_genes
plot(log(model1.counts[,1]), log(model1.counts[,3]))
parse_names
################EXPERIMENTAL DESIGN##################
experimental_design1 = data.frame(
  sample_names = col_names_counts,  # sample name
  treatment = factor(parse_names[,1]), # different tretments (3 in total)
  replicate = factor(parse_names[,2]),  # replicate1 or 2(each treatment has two replicate)
  barcode = factor(parse_names[,3]),   # 
  lane = factor(parse_names[,4])      # Whick lane on the Illumina flowcell.
)
write.csv(model1.counts, file="~/bio720/results/model1.counts.csv" )
write.csv(experimental_design1, file="~/bio720/results/experimental_design1.csv")

#######################################################################################################################
# Second MODEL : real samples and simulated samples (extended): test treatments
#######################################################################################################################
#ordering based on colnames

model2.counts <- data.frame(tot_count_sample, tot_count_simu1, tot_count_simu2, tot_count_simu3, tot_count_simu4, tot_count_simu5)
rownames(model2.counts) <- total_genes

parse_names2 <- strsplit(colnames(model2.counts), split="_")
parse_names2 <- matrix(unlist(parse_names2), nrow=36, ncol=3, byrow=T)

################EXPERIMENTAL DESIGN##################
experimental_design2 = data.frame(
  sample_names = colnames(model2.counts),  # sample name
  treatment = factor(parse_names2[,1]), # different tretments (3 in total)
  replicate = factor(parse_names2[,2]),  # replicate1 or 2(each treatment has two replicate)
  simulation = factor(parse_names2[,3])      # we had 5 round of simulation(S) and read samples(R)
)
write.csv(model2.counts, file="~/bio720/results/model2.counts.csv" )
write.csv(experimental_design2, file="~/bio720/results/experimental_design2.csv")



#######################################################################################################################


