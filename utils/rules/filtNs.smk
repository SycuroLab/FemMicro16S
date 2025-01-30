rule filterNsRaw:
    input:
        R1= config["input_dir"]+"/{sample}" + config["forward_read_suffix"] + config["compression_suffix"],
        R2= config["input_dir"]+"/{sample}" + config["reverse_read_suffix"] + config["compression_suffix"]
    params:
        dir=config["output_dir"]+"/filtN/"
    output:
        R1=config["output_dir"]+"/filtN/{sample}" + config["forward_read_suffix"] + config["compression_suffix"],
        R2=config["output_dir"]+"/filtN/{sample}" + config["reverse_read_suffix"] + config["compression_suffix"]
    conda:
        "dada2"
    script:
        "../scripts/dada2/filtN.R"

