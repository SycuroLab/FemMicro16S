rule rawReadCount:
    input:
        list_files["R1"].iloc[1].rsplit("/", 1)[0]
    output:
        config["output_dir"]+"/dada2/rawReadCount.txt"
    conda: "QC"
    params:
        format=config["compression_suffix"]
    shell:
        "seqkit stats -a {input}/*{params} > {output}"



rule numSeqParse:
    input:rules.rawReadCount.output
    output:config["output_dir"]+"/dada2/Parsed_rawReadCount.txt"
    params:
        R1= config["forward_read_suffix"]+config["compression_suffix"],
        R2= config["reverse_read_suffix"]+config["compression_suffix"]
    shell:
        """
        sed 's/,//g' {input} | sed -n '1p;0~2p' | awk '{{print $1,$4}}' | rev | cut -d'/' -f 1 | rev | sed 's/{params.R1}\\|{params.R2}$//'  | tr ' ' '\\t' > {output}
        """



#Used to combine read counts from each step of dada2
rule combineReadCounts:
    input:
        config["output_dir"]+"/dada2/Parsed_rawReadCount.txt",
        config["output_dir"]+"/dada2/Nreads_filtered.txt",
        config["output_dir"]+"/dada2/Nreads_with_chimeras.txt",
        config["output_dir"]+"/dada2/Nreads_nochimera.txt"
    output:
        config["output_dir"]+"/dada2/Nreads.tsv",
    run:
        import pandas as pd
        D= pd.read_table(input[0],index_col=0)
        D = D.join(pd.read_table(input[1], index_col=0), how='outer')
        D = D.join(pd.read_table(input[2], index_col=0), how='outer')
        D = D.join(pd.read_table(input[3], index_col=0), how='outer')
        D.to_csv(output[0],sep='\t')


