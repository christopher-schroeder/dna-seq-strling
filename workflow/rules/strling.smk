
# rule vembrane_filter:
#     input:
#         vcf="results/strling/vcf/{sample}.annotated.vcf"
#     output:
#         vcf="results/strling/vcf/{sample}.filtered.vcf"
#     conda:
#         "../envs/vembrane.yaml"
#     params:
#         expression="10**(-QUAL/10) < 0.1"
#     log:
#         "logs/vembrane-filter/{sample}.log"
#     shell:
#         "cat {input.vcf} | vembrane filter \"{params.expression}\" > {output.vcf} 2> {log}" #"1-10**(-INFO['PROB_DENOVO']/10)"

#

# rule snpeff:
#     input:
#         calls="results/strling/vcf/{group}/{group}.all.bcf",
#         db="resources/snpeff/GRCh38.86"
#     output:
#         calls="results/strling/vcf/{group}/{group}.all.annotated.bcf",
#         stats="results/strling/snpeff/{group}.html",
#         csvstats="results/strling/snpeff/{group}.csv"
#     log:
#         "logs/snpeff/{group}.log"
#     params:
#         extra="-nodownload"
#     resources:
#         mem_mb=4000
#     wrapper:
#         "0.68.0/bio/snpeff/annotate"


# rule snpeff_download:
#     output:
#         # wildcard {reference} may be anything listed in `snpeff databases`
#         directory("resources/snpeff/GRCh38.86")
#     log:
#         "logs/snpeff/download/GRCh38.86.log"
#     params:
#         reference="GRCh38.86"
#     resources:
#         mem_mb=4096
#     wrapper:
#         "0.68.0/bio/snpeff/download"


rule visualize:
    input:
        bcf="results/strling/vcf/{group}/{group}.all.annotated.bcf"
    output:
        "results/strling/plots/{group}.pdf"
    shell:
        "python workflow/scripts/visualize.py {input.bcf} {output}"


def merge_command(wc):
    if len(get_group_samples(wc.group)) == 1:
        return "cat"
    else:
        return "bcftools merge -m none -Ob"


rule annotate_quantile:
    input:
        genotypes=lambda wc: expand("results/strling/call/{{group}}/{sample}-genotype.txt", sample=set(samples.sample_name) - set(get_group_samples(wc.group))),
        bcf="results/strling/vcf/{group}/{group}.all.bcf"
    output:
        "results/strling/vcf/{group}/{group}.all.annotated.bcf"
    shell:
        "bcftools view {input.bcf} | python workflow/scripts/annotate_quantiles2.py {input.genotypes} | bcftools view -Ob > {output}"


rule merge_bcf:
    input:
        bcf=lambda wc: expand("results/strling/vcf/{{group}}/{sample}.bcf", sample=get_group_samples(wc.group)),
        index=lambda wc: expand("results/strling/vcf/{{group}}/{sample}.bcf.csi", sample=get_group_samples(wc.group))
    output:
        "results/strling/vcf/{group}/{group}.all.bcf"
    params:
        command=merge_command
    shell:
        "{params.command} {input.bcf} | bcftools norm -m -both > {output}"


rule index_bcf:
    input:
        "{x}.bcf"
    output:
        "{x}.bcf.csi"
    shell:
        "bcftools index {input}"


rule vcf_to_bcf:
    input:
        vcf="{x}.vcf"
    output:
        bcf="{x}.bcf"
    shell:
        "bcftools sort -Ob {input.vcf} | bcftools norm -Ob -m-any > {output.bcf}"


rule strling_to_vcf:
    input:
        tsv="results/strling/outlier/{group}/{sample}.STRs.tsv"
    output:
        vcf="results/strling/vcf/{group}/{sample}.vcf"
    params:
        samplename = "{sample}"
    script:
        "../scripts/strling_to_vcf.py"


# rule annotate_quantile:
#     input:
#         outliers="results/strling/outlier/{group}/{sample}.STRs.tsv",
#         quantiles="results/strling/quantiles/{group}.txt"
#     output:
#         "results/strling/outlier_quantiles/{group}/{sample}.STRs.tsv"
#     script:
#         "../scripts/annotate_quantiles.py"


rule strling_outlier:
    input:
        genotype=expand("results/strling/call/{{group}}/{sample}-genotype.txt", sample=samples.sample_name),
        unplaced=expand("results/strling/call/{{group}}/{sample}-unplaced.txt", sample=samples.sample_name)
    output:
        outliers=expand("results/strling/outlier/{{group}}/{sample}.STRs.tsv", sample=samples.sample_name),
        unplaced="results/strling/outlier/{group}/unplaced.tsv"
    params:
        prefix="results/strling/outlier/{group}/"
    shell:
        "strling-outliers.py --genotypes {input.genotype} --unplaced {input.unplaced} --out {params.prefix}"


# rule calculate_quantiles:
#     input:
#         genotypes=lambda wc: expand("results/strling/call/{{group}}/{sample}-genotype.txt", sample=set(samples.sample_name) - set(get_group_samples(wc.group)))
#     output:
#         quantiles="results/strling/quantiles/{group}.txt"
#     script:
#         "../scripts/get_quantiles.py"


rule strling_call:
    input:
        bam=get_bam,
        bai=get_bai,
        bin="results/strling/extract/{sample}.bin",
        reference="resources/genome.fasta",
        fai="resources/genome.fasta.fai",
        bounds="results/strling/merge/{group}-bounds.txt"
    output:
        bounds="results/strling/call/{group}/{sample}-bounds.txt",
        genotype="results/strling/call/{group}/{sample}-genotype.txt",
        unplaced="results/strling/call/{group}/{sample}-unplaced.txt",
    params:
        prefix="results/strling/call/{group}/{sample}"
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/call/{group}/{sample}.log"
    shell:
        "strling call -f {input.reference} -b {input.bounds} -o {params.prefix} {input.bam} {input.bin} 2> {log}"


rule strling_merge:
    input:
        bins=lambda w: expand("results/strling/extract/{sample}.bin", sample=get_group_samples(w.group)),
        reference="resources/genome.fasta",
        fai="resources/genome.fasta.fai",
    output:
        bounds="results/strling/merge/{group}-bounds.txt"
    params:
        prefix="results/strling/merge/{group}"
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/merge/{group}.log"
    shell:
        "strling merge -f {input.reference} -o {params.prefix} {input.bins} 2> {log}"


rule strling_extract:
    input:
        bam=get_bam,
        bai=get_bai,
        reference="resources/genome.fasta",
        fai="resources/genome.fasta.fai",
        index="resources/genome.fasta.str" # optional
    output:
        bin="results/strling/extract/{sample}.bin",
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/extract/{sample}.log"
    shell:
        "strling extract -f {input.reference} -g {input.index} {input.bam} {output.bin} 2> {log}"


rule strling_index:
    input:
        reference="resources/genome.fasta"
    output:
        "resources/genome.fasta.str",
        "resources/genome.fasta.fai",
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/index/genome.log"
    shell:
        "strling index {input.reference} 2> {log}"
