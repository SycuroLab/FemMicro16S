suppressMessages(library(dada2))


seqtab.all= readRDS(snakemake@input[['seqtab']]) # seqtab.all


# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab.all, method=snakemake@config[['chimera_method']], multithread=snakemake@threads)

A <- sum(seqtab) / sum(seqtab.all)

removed_percent <- (1 - A) * 100

message <- paste0(
  "Abundances of removed chimeric sequence variants was ",
  round(removed_percent, 2), "% of the merged sequence reads"
)
print(message)

saveRDS(seqtab, snakemake@output[['rds']])
write.csv(seqtab, snakemake@output[['csv']])


track <- rowSums(seqtab)
names(track) <- row.names(seqtab)

write.table(track,col.names = c("nonchim"),
            snakemake@output[['nreads']],sep='\t')
