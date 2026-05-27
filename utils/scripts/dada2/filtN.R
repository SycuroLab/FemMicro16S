suppressMessages(library(dada2))

Dir <- snakemake@params[["dir"]]

if (!dir.exists(Dir)) {
  dir.create(Dir, recursive = TRUE)
}

fnFs <- snakemake@input[["R1"]]
fnRs <- snakemake@input[["R2"]]

fnFs.filtN <- snakemake@output[["R1"]]
fnRs.filtN <- snakemake@output[["R2"]]

percent_phix_file <- snakemake@output[["percent_phix"]]

tmpF <- paste0(fnFs.filtN, ".no_phix_check.tmp.gz")
tmpR <- paste0(fnRs.filtN, ".no_phix_check.tmp.gz")

out_no_phix_removal <- filterAndTrim(
  fnFs, tmpF,
  fnRs, tmpR,
  maxN = 0,
  rm.phix = FALSE,
  multithread = TRUE
)

out_with_phix_removal <- filterAndTrim(
  fnFs, fnFs.filtN,
  fnRs, fnRs.filtN,
  maxN = 0,
  rm.phix = TRUE,
  multithread = TRUE
)

percent_phix <- data.frame(
  sample = basename(fnFs),
  reads_input = out_no_phix_removal[, 1],
  reads_after_maxN_no_phix_removal = out_no_phix_removal[, 2],
  reads_after_maxN_with_phix_removal = out_with_phix_removal[, 2],
  percent_phix = 100 * (
    out_no_phix_removal[, 2] - out_with_phix_removal[, 2]
  ) / out_no_phix_removal[, 1]
)

write.table(
  percent_phix,
  percent_phix_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

file.remove(tmpF)
file.remove(tmpR)

