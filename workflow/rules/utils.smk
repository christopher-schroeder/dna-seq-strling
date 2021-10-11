rule download_genehancer:
    input:
        fai="resources/genome.fasta.fai"
    output:
        "results/resources/genehancer.tsv"
    params:
        hg="hg38"
    log:
        "logs/download_genehancer.log"
    cache: True
    script:
        "../scripts/dl_genehancer.py"
