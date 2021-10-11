rule download_genehancer:
    input:
        fai="results/resources/genome.fasta.fai"
    output:
        "results/resources/genehancer.tsv"
    params:
        hg="hg38"
    log:
        "logs/download_genehancer.log"
    cache: True
    script:
        "../scripts/dl_genehancer.py"


rule index_bcf:
    input:
        "{x}.bcf"
    output:
        "{x}.bcf.csi"
    shell:
        "bcftools index {input}"