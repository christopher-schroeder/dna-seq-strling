import pandas as pd
from snakemake.utils import min_version

min_version("5.12.0")

configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/expansionhunter.smk"

# targets = []
# for group in groups_no_control:
#     targets.extend(expand("results/strling/vcf/{group}/{group}.all.annotated.bcf", group=[group]))

rule all:
    input:
        expand("results/eh/calls/{sample}.json", sample=samples["sample_name"])