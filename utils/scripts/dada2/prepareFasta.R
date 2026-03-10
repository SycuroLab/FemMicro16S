suppressMessages(library("dplyr"))

# Create a function to write FASTA entries
write_fasta <- function(df, file) {
  sink(file)
  for (i in 1:nrow(df)) {
    cat(paste0(">", df[i,1], "\n"))
    cat(paste0(df[i,2], "\n"))
  }
  sink()
}


df<-read.csv(snakemake@input[['file']], sep = ",", header = TRUE)


df_selected <- df %>%
  select(asv_num, asv_seq, matches("GTDB", ignore.case = TRUE))



# Define the function to format taxonomy columns
format_taxonomy <- function(row) {
  paste0(
    "k__", row[1], ";",
    "p__", row[2], ";",
    "c__", row[3], ";",
    "o__", row[4], ";",
    "f__", row[5], ";",
    "g__", row[6], ";",
    "s__", row[7]
  )
}

# Apply the function within the mutate step
df_selected <- df %>%
  select(asv_num, asv_seq, matches("GTDB", ignore.case = TRUE)) %>%
  mutate(taxonomy = paste0(asv_num,";",apply(.[, 3:9], 1, format_taxonomy)))


# Write the data to a FASTA file with asv ids
write_fasta(df_selected %>% select("asv_num","asv_seq"), snakemake@output[['id_fas']])


# Write the data to a FASTA file with asv ids and taxonomy
write_fasta(df_selected %>% select("taxonomy","asv_seq"), snakemake@output[['id_tax']])
