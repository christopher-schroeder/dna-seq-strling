import pandas as pd
from snakemake.utils import min_version
min_version("5.12.0")

configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/strling.smk"
include: "rules/table.smk"
include: "rules/utils.smk"

# targets = []
# for group in groups_no_control:
#     targets.extend(expand("results/strling/vcf/{group}/{group}.all.annotated.bcf", group=[group]))


rule all:
    input:
        expand("results/strling/plots/{group}.pdf", group=groups_no_control),
        "results/resources/genehancer.tsv"
        #expand("results/strling/tables/{group}/{group}.xlsx", group=groups_no_control)


rule only_bai_linking:
    input:
        expand("{x}.bai", x=samples["bam"])
