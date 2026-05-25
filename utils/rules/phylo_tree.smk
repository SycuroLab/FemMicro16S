rule multipleAlign:
    input:
        seqtab =rules.removeChimeras.output.rds
    output:
        seqfasta=config["output_dir"]+"/fasta_files/ASVs_seqs.fasta",
        alignment=config["output_dir"]+"/phylogeny/ASV_aligned.fasta"
    threads:
        config['threads']
    singularity:
        "apptainer/dada2-1.0.0.sif"
    script:
        "../scripts/dada2/alignment.R"

   

rule newickTree:
    input:
         rules.multipleAlign.output.alignment
    output:
        config["output_dir"]+"/phylogeny/ASV_tree.nwk"
    threads:
        config['threads']
    singularity:
        "apptainer/fastree_mafft-1.0.0.sif"
    shell:
        """
        export OMP_NUM_THREADS={threads} && FastTreeMP -nt -gamma -spr 4 {input} > {output} || fasttree -nt -gamma -spr 4  {input} > {output}
        """
