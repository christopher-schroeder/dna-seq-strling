
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

rule visualize:
    input:
        bcf="results/strling/vcf/{experiment}/{experiment}.all.annotated.bcf"
    output:
        report(
            directory("results/strling/plots/{experiment}"),
            patterns=["{chrom}-{left}-{right}-{motif}-{gene}.pdf"],
            caption="../report/plots.rst",
            category="STR plots",
            subcategory="{experiment}",
        ),
    conda:
        "../envs/visualize.yaml"
    shell:
        "python workflow/scripts/visualize_single.py {input.bcf} {output}"


def merge_command(wc):
    if len(get_experiment_samples(wc.experiment)) == 1:
        return "cat"
    else:
        return "bcftools merge -m none -Ob"


rule annotate_quantile:
    input:
        genotypes=lambda wc: expand("results/strling/call/{{experiment}}/{sample}-genotype.txt", sample=get_experiment_samples(wc.experiment, case=False, control=True)),
        bcf="results/strling/vcf/{experiment}/{experiment}.all.bcf"
    output:
        "results/strling/vcf/{experiment}/{experiment}.all.q.bcf"
    shell:
        "bcftools view {input.bcf} | python workflow/scripts/annotate_quantiles2.py {input.genotypes} | bcftools view -Ob > {output}"


rule merge_bcf:
    input:
        bcf=lambda wc: expand("results/strling/vcf/{{experiment}}/{sample}.bcf", sample=get_experiment_samples(wc.experiment, case=True, control=False)),
        index=lambda wc: expand("results/strling/vcf/{{experiment}}/{sample}.bcf.csi", sample=get_experiment_samples(wc.experiment, case=True, control=False))
    output:
        "results/strling/vcf/{experiment}/{experiment}.all.bcf"
    params:
        command=merge_command
    shell:
        "{params.command} {input.bcf} | bcftools norm -m -both > {output}"


rule vcf_to_bcf:
    input:
        vcf="results/strling/vcf/{experiment}/{sample}.vcf"
    output:
        bcf="results/strling/vcf/{experiment}/{sample}.bcf"
    shell:
        "bcftools sort -Ob {input.vcf} | bcftools norm -Ob -m-any > {output.bcf}"


rule strling_to_vcf:
    input:
        tsv="results/strling/outlier/{experiment}/"
    output:
        vcf="results/strling/vcf/{experiment}/{sample}.vcf"
    params:
        tsv = "results/strling/outlier/{experiment}/{sample}.STRs.tsv",
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
        genotype=lambda wc: expand("results/strling/call/{{experiment}}/{sample}-genotype.txt", sample=get_experiment_samples(wc.experiment, case=True, control=True)),
        unplaced=lambda wc: expand("results/strling/call/{{experiment}}/{sample}-unplaced.txt", sample=get_experiment_samples(wc.experiment, case=True, control=True)),
    output:
        outliers=directory("results/strling/outlier/{experiment}"),
        unplaced="results/strling/outlier/{experiment}/unplaced.tsv"
    conda:
        "../envs/strling.yaml"
    params:
        prefix="results/strling/outlier/{experiment}/"
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
        reference=genome,
        fai=genome + ".fai",
        bounds="results/strling/merge/{experiment}-bounds.txt"
    output:
        bounds="results/strling/call/{experiment}/{sample}-bounds.txt",
        genotype="results/strling/call/{experiment}/{sample}-genotype.txt",
        unplaced="results/strling/call/{experiment}/{sample}-unplaced.txt",
    params:
        prefix="results/strling/call/{experiment}/{sample}"
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/call/{experiment}/{sample}.log"
    shell:
        "/vol/nano/christo/STRling/strling call -f {input.reference} -b {input.bounds} -o {params.prefix} {input.bam} {input.bin} 2> {log}"


rule strling_merge:
    input:
        bins=lambda w: expand("results/strling/extract/{sample}.bin", sample=get_experiment_samples(w.experiment, case=True, control=False)),
        reference=genome,
        fai=genome + ".fai",
    output:
        bounds="results/strling/merge/{experiment}-bounds.txt"
    params:
        prefix="results/strling/merge/{experiment}"
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/merge/{experiment}.log"
    shell:
        "/vol/nano/christo/STRling/strling merge -f {input.reference} -o {params.prefix} {input.bins} 2> {log}"


rule strling_extract:
    input:
        bam=get_bam,
        bai=get_bai,
        reference=genome,
        fai=genome + ".fai",
        index=genome + ".str"
    output:
        bin="results/strling/extract/{sample}.bin",
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/extract/{sample}.log"
    shell:
        "/vol/nano/christo/STRling/strling extract -f {input.reference} -g {input.index} {input.bam} {output.bin} 2> {log}"


rule strling_index:
    input:
        reference=genome
    output:
        str_index=genome + ".str",
        # genome + ".fai",
    conda:
        "../envs/strling.yaml"
    log:
        "logs/strling/index/genome.log"
    cache: True
    shell:
        "/vol/nano/christo/STRling/strling index {input.reference} -g {output.str_index} 2> {log}"


ruleorder: strling_index > genome_faidx