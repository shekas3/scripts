setwd("c:/Users/guest6/Desktop/FinalProject/counts_htseq/")
path = "c:/Users/guest6/Desktop/FinalProject/counts_htseq/"
file.names <- dir(path, pattern = ".txt")
file.counts <- list ()
All.samples <- list ()
for(i in 1:length(file.names)){
  file.counts[[ file.names[i] ]] <- read.table(file.names[i], header = FALSE, sep = "\t")
  All.samples <- data.frame(file.counts)
  }
comp.samples <- All.samples[,-c(3,5,7,9,11)]
cols1 <- c("genes", "SD", "SDR", "SRW.1", "SRW.2", "SWW.3", "SWW")
colnames(comp.samples) <- cols1
#################################################
plot(log(comp.samples$SWW),log(comp.samples$SWW.3) )
plot(log(comp.samples$SD), log(comp.samples$SDR))
plot(log(comp.samples$SRW.1), log(comp.samples$SRW.2))
#######################################################
#Total number of reads that mapped to a gene or total number of counts for
#sample: this is necesarrly for simulation
sapply(comp.samples, function(x) sum(as.numeric(x)) )

library(splitstackshape)
total.counts.SDR <- expandRows(comp.samples, "SDR")

