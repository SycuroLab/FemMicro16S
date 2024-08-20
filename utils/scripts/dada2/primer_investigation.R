suppressMessages(library(dada2))
suppressMessages(library(Biostrings))
suppressMessages(library(ShortRead))


Dir <- snakemake@params[["dir"]]

if (!dir.exists(Dir)) {
  dir.create(Dir)
}



####Checking presence of primers before removing the primers

fnFs.filtN <- snakemake@input[['R1']]
fnRs.filtN <- snakemake@input[['R2']]

fnFs.cut <- snakemake@input[['cut1']]
fnRs.cut <- snakemake@input[['cut2']]


FWD <- snakemake@config[["fwd_primer"]]
REV <- snakemake@config[["rev_primer"]]


fwd_pattern <- paste0(snakemake@config[["forward_read_suffix"]],snakemake@config[["compression_suffix"]])
rev_pattern <- paste0(snakemake@config[["reverse_read_suffix"]],snakemake@config[["compression_suffix"]])


allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)




primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}



primer_status_bf <- NULL

for (i in 1:length(fnFs.filtN)) {
  temp <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[i]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[i]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[i]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[i]])
  )
  sample <- gsub(x = basename(fnFs.filtN[i]),pattern = paste0(fwd_pattern,"|",rev_pattern),replacement = "")
  row <- cbind(sample, temp)
  primer_status_bf <- rbind(primer_status_bf, row)
}


write.csv(primer_status_bf,snakemake@output[["primer_status_bf"]])


####Checking presence of primers after removing the primers

primer_status_af=NULL
	
for (i in 1:length(fnFs.cut)) {
  temp<-rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[i]]), 
              FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[i]]), 
              REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[i]]), 
              REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[i]]))
  sample <- gsub(x = basename(fnFs.cut[i]),pattern = paste0(fwd_pattern,"|",rev_pattern),replacement = "") 
  row=cbind(sample,temp)
  primer_status_af<-rbind(primer_status_af,row)
}



write.csv(primer_status_af,snakemake@output[["primer_status_af"]])
