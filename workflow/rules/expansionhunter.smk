rule expansionhunter:
    input:
        bam=get_bam,
        bai=get_bai,
        reference="resources/genome.fasta",
        catalogue="results/eh/catalogue/catalogue.json"
    output:
        json="results/eh/calls/{sample}.json",
        vcf="results/eh/calls/{sample}.vcf",
        bamlet="results/eh/calls/{sample}_realigned.bam"
    params:
        prefix="results/eh/calls/{sample}"
    log:
        "logs/eh/{sample}.log"
    conda:
        "../envs/expansionhunter.yaml"
    shell:
        """
        (ExpansionHunter --reads {input.bam} \
            --reference {input.reference} \
            --variant-catalog {input.catalogue} \
            --output-prefix {params.prefix}) > {log} 2>&1
        """


rule create_catalogue:
    input:
        regions="config/regions.tsv"
    output:
        catalogue="results/eh/catalogue/catalogue.json"
    script:
        "../scripts/create_catalogue.py"