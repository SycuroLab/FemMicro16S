rule vsearchGTDB:
    input:
        fas=rules.fasta_file.output.seqfasta
    output:
        output1= config["output_dir"]+"/vsearch/GTDB/Vsearch_GTDB_selected.tsv",
        output2= config["output_dir"]+"/vsearch/GTDB/Vsearch_GTDB_raw.tsv"
    threads:
        config['threads']
    conda:
        "vsearch"
    log:
        config["output_dir"]+"/vsearch/logs/GTDB_logfile.txt"
    params:
        db=config["vsearch_DBs"]["GTDB"],
        id=config["Identity"],
        maxaccepts=config["Maxaccepts"]
    shell: 
        "time vsearch --usearch_global {input.fas} --db {params.db} --notrunclabels --userout {output.output1} --userfields query+target+id --uc {output.output2} --id {params.id} --iddef 0 --log {log} --threads {threads} --uc_allhits --maxaccepts {params.maxaccepts} --top_hits_only --strand both --gapopen '*' "



rule separate_vsearch_hits:
    input:
        GTDB_file=rules.vsearchGTDB.output.output1,
        annotation=rules.combining_annotations.output.table
    output:
        All_ASVs=config["output_dir"]+"/vsearch/GTDB/All_results_GTDB.tsv",
        GTDB_noHIT=config["output_dir"]+"/vsearch/GTDB/no_hits_GTDB.fasta"
    conda:
        "dada2"  
    script:
        "../scripts/vsearch/vsearch.R"



rule vsearchURE:
    input:
        fas=rules.separate_vsearch_hits.output.GTDB_noHIT
    output:
        output1= config["output_dir"]+"/vsearch/URE/Vsearch_URE_selected.tsv",
        output2= config["output_dir"]+"/vsearch/URE/Vsearch_URE_raw.tsv"
    threads:
        config['threads']
    conda:
        "vsearch"
    log:
        config["output_dir"]+"/vsearch/logs/URE_logfile.txt"
    params:
        db=config["vsearch_DBs"]["URE"],
        id=config["Identity"],
        maxaccepts=config["Maxaccepts"]
    shell:
        """
        if [[ "{config[URE_after_GTDB]}" == True ]]; then
            time vsearch --usearch_global {input.fas} --db {params.db}   --notrunclabels --userout {output.output1} \
            --userfields query+target+id --uc {output.output2} --id {params.id} --iddef 0 --log {log} \
            --threads {threads} --uc_allhits --maxaccepts {params.maxaccepts} --top_hits_only --strand both --gapopen '*'
        else
            echo "Rule 'vsearchURE' is not executed."
        fi
        """



rule vsearchParse:
    input:
        GTDB_df=rules.separate_vsearch_hits.output.All_ASVs,
        annotation=rules.combining_annotations.output.table,
        URE_df=lambda wildcards: [rules.vsearchURE.output.output1] if config.get("URE_after_GTDB", True) else []
    output:
        SemiParsed_uncollapsed=config["output_dir"]+"/vsearch/Final_uncollapsed_output.tsv",
        parsed_collapsed_GTDB_URE=config["output_dir"]+"/vsearch/Final_colapsed_output.tsv",
        Vsearch_final=config["output_dir"]+"/taxonomy/vsearch_tables/Vsearch_output.tsv",
        merged_final=config["output_dir"]+"/taxonomy/final_merged_tables/vsearch_dada2_merged.tsv"
    threads:
        config['threads']
    conda:
        "dada2"   
    script:
        "../scripts/vsearch/vsearch_parsing.R"
