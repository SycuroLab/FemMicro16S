rule combining_annotations:
    input:
        expand(config["output_dir"] + "/taxonomy/dada2_tables/{ref}_RDP.tsv", ref=config['RDP_dbs'].keys()),
        seqs=rules.removeChimeras.output.rds
    output:
        table=config["output_dir"]+"/taxonomy/dada2_tables/"+"dada2_all_databases_merged.csv"
    conda:
        "dada2"
    script:
        "../scripts/dada2/annotation_output_dada2.R"



rule prepareFasta:
    input:
        file=rules.combining_annotations.output.table
    output:
        id_fas=config["output_dir"]+"/fasta_files/ASVs_id.fasta",
        id_tax=config["output_dir"]+"/fasta_files/ASVs_tax.fasta"
    conda:
        "dada2"
    script:
        "../scripts/dada2/prepareFasta.R"
