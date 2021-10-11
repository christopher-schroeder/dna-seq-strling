rule get_genome:
    output:
        "results/resources/genome.fasta"
    log:
        "logs/get-genome.log"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        "results/resources/genome.fasta"
    output:
        "results/resources/genome.fasta.fai"
    log:
        "logs/genome-faidx.log"
    cache: True
    wrapper:
        "0.59.2/bio/samtools/faidx"


rule genome_dict:
    input:
        "results/resources/genome.fasta"
    output:
        "results/resources/genome.dict"
    log:
        "logs/samtools/create_dict.log"
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "
