#Required python modules
import os
import pandas as pd


##Input file
configfile: "config.yaml"
list_files = pd.read_table(config['list_files'],index_col=0)
SAMPLES = list(list_files.index)




myoutput = list()

if config['primer_removal'] == True and config['primer_investigation'] == True:
    myoutput.append(config["output_dir"] + "/seqkit_samples/" + "temp_primerRMV.txt")
    myoutput.append(config["output_dir"] + "/primer_status/primer_existance_raw.csv")
    myoutput.append(config["output_dir"] + "/primer_status/primer_existance_trimmed.csv")



rule all:
    input:
        config["output_dir"]+"/figures/length_distribution/Sequence_Length_distribution.png",
        config["output_dir"]+"/figures/quality/afterQCQualityPlots"+ config["forward_read_suffix"]+".png",
        config["output_dir"]+"/figures/quality/afterQCQualityPlots"+ config["reverse_read_suffix"]+".png",
        config["output_dir"]+"/figures/quality/rawFilterQualityPlots"+ config["forward_read_suffix"]+".png",
        config["output_dir"]+"/figures/quality/rawFilterQualityPlots"+ config["reverse_read_suffix"]+".png",
        config["output_dir"]+"/figures/quality/afterdada2FilterQualityPlots"+ config["forward_read_suffix"]+".png",
        config["output_dir"]+"/figures/quality/afterdada2FilterQualityPlots"+ config["reverse_read_suffix"]+".png",
        expand(config["output_dir"]+"/taxonomy/dada2_tables/{ref}_RDP.tsv",ref= config['RDP_dbs'].keys()),
        config["output_dir"]+"/dada2/Nreads.tsv",
        config["output_dir"]+"/dada2/Nreads_filtered.txt",
        config["output_dir"]+"/dada2/percent_phix.txt",
        config["output_dir"]+"/multiqc_filt/multiqc_report_filtered.html",
        config["output_dir"]+"/multiqc_raw/multiqc_report_raw.html",
        config["output_dir"]+"/seqkit_samples/"+"temp_raw.txt",
        config["output_dir"]+"/seqkit_samples/"+"temp_dada2.txt",
        config["output_dir"]+"/seqkit_samples/"+"temp_cutadapt.txt",
        config["output_dir"]+"/QC_html_report/"+"final_qc_report.html",
        config["output_dir"]+"/taxonomy/dada2_tables/"+"dada2_all_databases_merged.csv",
        config["output_dir"]+"/vsearch/Final_uncollapsed_output.tsv",
        config["output_dir"]+"/vsearch/Final_colapsed_output.tsv",
        config["output_dir"]+"/taxonomy/vsearch_tables/Vsearch_output.tsv",
        config["output_dir"]+"/taxonomy/final_merged_tables/vsearch_dada2_merged.tsv",
        config["output_dir"]+"/fasta_files/ASVs_id.fasta",
        config["output_dir"]+"/fasta_files/ASVs_tax.fasta",
        config["output_dir"]+"/fasta_files/ASVs_seqs.fasta",
        myoutput


##path to where different snakemake rule files are saved

include: "utils/rules/filtNs.smk"
include: "utils/rules/qc_cutadapt.smk"
include: "utils/rules/dada2.smk"
include: "utils/rules/fasta_generation.smk"
include: "utils/rules/readCount.smk"
include: "utils/rules/seqkit_length_report.smk"
include: "utils/rules/annotation_output_dada2.smk"
include: "utils/rules/vsearch.smk"
include: "utils/rules/qc_report.smk"
