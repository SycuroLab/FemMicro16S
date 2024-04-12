rule combining_annotations:
    input:
        config["output_dir"]+"/taxonomy/dada2_tables/URE_RDP.tsv",
        config["output_dir"]+"/taxonomy/dada2_tables/GTDB_RDP.tsv",
        config["output_dir"]+"/taxonomy/dada2_tables/RDP_RDP.tsv",
        config["output_dir"]+"/taxonomy/dada2_tables/Silva_RDP.tsv",
        seqs=rules.removeChimeras.output.rds 
    output:
        table=config["output_dir"]+"/taxonomy/dada2_tables/"+"dada2_all_databases_merged.tsv"
    conda:
        "dada2"
    script:
        "../scripts/dada2/annotation_output_dada2.R"
