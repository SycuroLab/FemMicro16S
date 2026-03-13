suppressMessages(library(dada2))


fnFs <- snakemake@input[['R1']]
fnRs <- snakemake@input[['R2']]

filtFs <- snakemake@output[['R1']]
filtRs <- snakemake@output[['R2']]


track.filt <- filterAndTrim(fnFs,filtFs, 
                            fnRs,filtRs, 
#                            truncLen= snakemake@config[["truncLen"]],
                            maxN=0,
                            maxEE=snakemake@config[["maxEE"]], 
                            truncQ=snakemake@config[["truncQ"]],
                            compress=TRUE,
                            verbose=TRUE,
                            multithread=snakemake@threads,
                            rm.phix=TRUE)

row.names(track.filt) <- snakemake@params[["samples"]]

colnames(track.filt) = c('afterCutadapt','filtered')

track.filt<-data.frame(track.filt)
track.filt <- track.filt[track.filt$filtered > 0, ]

write.table(track.filt,snakemake@params[['nread']],  sep='\t')
