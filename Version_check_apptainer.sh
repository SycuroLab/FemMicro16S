#!/usr/bin/env bash
# Version_check_apptainer_simple.sh

set -u
set -o pipefail

echo "========== Host tools =========="
echo "Conda:      $(conda --version 2>&1 | head -n 1)"
echo "Mamba:      $(mamba --version 2>&1 | head -n 1)"
echo "Pip:        $(pip --version 2>&1 | head -n 1)"

echo
echo "========== snakemake-1.0.0.sif =========="
echo "Snakemake (snakemake-1.0.0.sif): $(apptainer exec apptainer/snakemake-1.0.0.sif snakemake --version 2>&1 | head -n 1)"
echo "PyYAML (snakemake-1.0.0.sif):    $(apptainer exec apptainer/snakemake-1.0.0.sif python -c 'import yaml, sys; sys.stdout.write(getattr(yaml, "__version__", "?"))' 2>&1 | head -n 1)"

echo
echo "========== dada2-1.0.0.sif =========="
echo "R (dada2-1.0.0.sif):             $(apptainer exec apptainer/dada2-1.0.0.sif R --version 2>&1 | head -n 1)"

echo "dada2 (dada2-1.0.0.sif):         $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("dada2")),"\n"))' 2>&1 | tail -n 1)"
echo "DECIPHER (dada2-1.0.0.sif):      $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("DECIPHER")),"\n"))' 2>&1 | tail -n 1)"
echo "Biostrings (dada2-1.0.0.sif):    $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("Biostrings")),"\n"))' 2>&1 | tail -n 1)"
echo "limma (dada2-1.0.0.sif):         $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("limma")),"\n"))' 2>&1 | tail -n 1)"
echo "dplyr (dada2-1.0.0.sif):         $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("dplyr")),"\n"))' 2>&1 | tail -n 1)"
echo "ggplot2 (dada2-1.0.0.sif):       $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("ggplot2")),"\n"))' 2>&1 | tail -n 1)"
echo "gridExtra (dada2-1.0.0.sif):     $(apptainer exec apptainer/dada2-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("gridExtra")),"\n"))' 2>&1 | tail -n 1)"

# seqtk prints usage and exits non-zero; ignore failure
echo "seqtk (dada2-1.0.0.sif):         $(apptainer exec apptainer/dada2-1.0.0.sif seqtk 2>&1 | grep -i Version)"

echo
echo "========== qc-1.0.0.sif =========="
echo "Python (qc-1.0.0.sif):           $(apptainer exec apptainer/qc-1.0.0.sif python --version 2>&1 | head -n 1)"
echo "FastQC (qc-1.0.0.sif):           $(apptainer exec apptainer/qc-1.0.0.sif fastqc --version 2>&1 | head -n 1)"
echo "MultiQC (qc-1.0.0.sif):          $(apptainer exec apptainer/qc-1.0.0.sif multiqc --version 2>&1 | head -n 1)"

echo "pandas (qc-1.0.0.sif):           $(apptainer exec apptainer/qc-1.0.0.sif python3 -m pip show pandas | grep -i Version)"
echo "cutadapt (qc-1.0.0.sif):         $(apptainer exec apptainer/qc-1.0.0.sif cutadapt --version 2>&1 | head -n 1)"
echo "seqkit (qc-1.0.0.sif):           $(apptainer exec apptainer/qc-1.0.0.sif seqkit 2>&1 | grep -i 'Version:' | head -n 1 || true)"

echo
echo "========== fastree_mafft-1.0.0.sif =========="
# FastTree prints help + non-zero; ignore failure
echo "FastTree (fastree_mafft-1.0.0.sif): $(apptainer exec apptainer/fastree_mafft-1.0.0.sif FastTree -help 2>&1 | head -n 1 || true)"
echo "MAFFT (fastree_mafft-1.0.0.sif):      $(apptainer exec apptainer/fastree_mafft-1.0.0.sif mafft --version 2>&1 | head -n 1)"

echo
echo "========== rmd-1.0.0.sif =========="
# KronaTools version number between "KronaTools " and " - ktImportText"
echo "Krona ktImportText (rmd-1.0.0.sif): $(apptainer exec apptainer/rmd-1.0.0.sif ktImportText 2>&1 | grep -Po '(?<=KronaTools )[0-9.]+(?= - ktImportText)' | head -n 1 || true)"
echo "git (rmd-1.0.0.sif):                $(apptainer exec apptainer/rmd-1.0.0.sif git --version 2>&1 | head -n 1)"
echo "perl (rmd-1.0.0.sif):               $(apptainer exec apptainer/rmd-1.0.0.sif perl -e 'print $^V, "\n"' 2>&1 | head -n 1)"
echo "tar (rmd-1.0.0.sif):                $(apptainer exec apptainer/rmd-1.0.0.sif tar --version 2>&1 | head -n 1)"
echo "pandoc (rmd-1.0.0.sif):             $(apptainer exec apptainer/rmd-1.0.0.sif pandoc --version 2>&1 | head -n 1)"

echo "tidyverse (rmd-1.0.0.sif):          $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("tidyverse")),"\n"))' 2>&1 | tail -n 1)"
echo "kableExtra (rmd-1.0.0.sif):         $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("kableExtra")),"\n"))' 2>&1 | tail -n 1)"
echo "ggpubr (rmd-1.0.0.sif):             $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("ggpubr")),"\n"))' 2>&1 | tail -n 1)"
echo "DT (rmd-1.0.0.sif):                 $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("DT")),"\n"))' 2>&1 | tail -n 1)"
echo "RColorBrewer (rmd-1.0.0.sif):       $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("RColorBrewer")),"\n"))' 2>&1 | tail -n 1)"
echo "waterfalls (rmd-1.0.0.sif):         $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("waterfalls")),"\n"))' 2>&1 | tail -n 1)"
echo "plotly (rmd-1.0.0.sif):             $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("plotly")),"\n"))' 2>&1 | tail -n 1)"
echo "remotes (rmd-1.0.0.sif):            $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("remotes")),"\n"))' 2>&1 | tail -n 1)"
echo "phyloseq (rmd-1.0.0.sif):           $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("phyloseq")),"\n"))' 2>&1 | tail -n 1)"
echo "Biobase (rmd-1.0.0.sif):            $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("Biobase")),"\n"))' 2>&1 | tail -n 1)"
echo "Biostrings (rmd-1.0.0.sif):         $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("Biostrings")),"\n"))' 2>&1 | tail -n 1)"
echo "limma (rmd-1.0.0.sif):              $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages(cat(as.character(utils::packageVersion("limma")),"\n"))' 2>&1 | tail -n 1)"
echo "psadd (rmd-1.0.0.sif):              $(apptainer exec apptainer/rmd-1.0.0.sif Rscript --vanilla -e 'suppressPackageStartupMessages({if (requireNamespace("psadd", quietly=TRUE)) cat(as.character(utils::packageVersion("psadd")),"\n") 
else cat("NOT INSTALLED\n")})' 2>&1 | tail -n 1)"

echo
echo "========== vsearch-1.0.0.sif =========="
# vsearch prints banner + exits non-zero without args; ignore failure
echo "vsearch (vsearch-1.0.0.sif):       $(apptainer exec apptainer/vsearch-1.0.0.sif vsearch 2>&1 | head -n 1 || true)"

