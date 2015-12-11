setwd("~/bio720/counts_htseq/")
path = "~/bio720/counts_htseq/"
file.names <- dir(path, pattern = ".txt")
file.counts <- list ()
All.samples <- list ()
for(i in 1:length(file.names)){
  file.counts[[ file.names[i] ]] <- read.table(file.names[i], header = FALSE, sep = "\t")
  All.samples <- data.frame(file.counts)
}
tail(All.samples)
comp.samples <- All.samples[,-c(3,5,7,9,11)]
cols1 <- c("genes", "SDR", "SD", "SRW.1", "SRW.2", "SWW.3", "SWW")
colnames(comp.samples) <- cols1
tail(comp.samples)
head(comp.samples)
######################################################
plot(log(comp.samples$SWW),log(comp.samples$SWW.3) )
plot(log(comp.samples$SD), log(comp.samples$SDR))
plot(log(comp.samples$SRW.1), log(comp.samples$SRW.2))
#######################################################
#Total number of reads that mapped to a gene or total number of counts for
#sample: this is necesarrly for simulation
sapply(comp.samples, function(x) sum(as.numeric(x)) )

library(splitstackshape)
SDR.counts <- comp.samples[,c("genes", "SDR")]
tail(SDR.counts)
total.counts.SDR <- expandRows(SDR.counts, "SDR")
tail(total.counts.SDR)
resample.total.counts.SDR <- data.frame(total.counts.SDR[sample(150352759, replace = TRUE), ])
colnames(resample.total.counts.SDR) <- "genes"
colnames(resample.total.counts.SDR)
SDR.replicate3 <- data.frame(table(resample.total.counts.SDR))

head(SDR.simu)
head(SDR.counts)
tail(SDR.simu)
tail(SDR.counts)
tail(All.samples)
