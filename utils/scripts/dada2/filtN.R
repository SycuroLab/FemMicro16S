suppressMessages(library(dada2))

Dir <- snakemake@params[["dir"]]

if (!dir.exists(Dir)) {
  dir.create(Dir)
}

####Checking presence of primers before removing the primers

fnFs <- snakemake@input[['R1']]
fnRs <- snakemake@input[['R2']]

fnFs.filtN <- snakemake@output[['R1']] # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- snakemake@output[['R2']]

out<-filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, rm.phix=FALSE, multithread = TRUE) #This is to calculate phix reads percentage

out1<-filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, rm.phix=TRUE, multithread = TRUE) #These results are going to be used for downstream analysis


percent_phix=data.frame(percent_phix=100*(out[,2]-out1[,2])/out[,1])

write.table(percent_phix,snakemake@params[['percent_phix']],  sep='\t')

