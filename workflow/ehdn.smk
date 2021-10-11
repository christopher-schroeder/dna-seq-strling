import pandas as pd
from snakemake.utils import min_version
min_version("5.12.0")

configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/expansionhunterdenovo.smk"

rule all:
    input:
        expand("results/ehdn/outlier/{sample}.outlier_locus.tsv", sample=samples.sample_name)