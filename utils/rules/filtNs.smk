rule filterNsRaw:
    input:
        R1= expand(config["input_dir"]+"/{sample}" + config["forward_read_suffix"] + config["compression_suffix"],sample=SAMPLES),
        R2= expand(config["input_dir"]+"/{sample}" + config["reverse_read_suffix"] + config["compression_suffix"],sample=SAMPLES)
    params:
        dir=config["output_dir"]+"/filtN/",
        percent_phix= config["output_dir"]+"/filtN/percent_phix.txt"
    output:
        R1=config["output_dir"]+"/filtN/{sample}" + config["forward_read_suffix"] + config["compression_suffix"],
        R2=config["output_dir"]+"/filtN/{sample}" + config["reverse_read_suffix"] + config["compression_suffix"]
    singularity:
        "apptainer/dada2-1.0.0.sif"
    script:
        "../scripts/dada2/filtN.R"

