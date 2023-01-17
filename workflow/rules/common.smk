import glob
import pandas as pd
import os
import yaml

samples = pd.read_csv(config["samples"], sep="\t", comment='#', dtype={"sample_name": str, "bam": str, "group":str}).set_index("sample_name", drop=False).sort_index()
groups = samples["group"].unique()

groups_no_control = groups[groups != "control"]
samples_no_control = samples[samples["group"] != "control"]
experiments = dict()

for filename in glob.glob("config/experiments/*.yaml"):
    name = os.path.basename(filename)[:-len(".yaml")]
    with open(filename, "r") as stream:
        experiments[name] = yaml.safe_load(stream)

wildcard_constraints:
    group="|".join(samples["group"].unique()),
    sample="|".join(samples["sample_name"]),
    experiment="|".join(experiments)

datatype = "dna"
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
genome = f"resources/genome.{datatype}.{species}.{build}.{release}.fasta"

def get_annotations_extra(wildcards, input):
    if annotations:
        custom_str = " ".join(f"--custom {path},{prefix},vcf,exact,,{','.join(ann['fields'])}" for (prefix, ann), path in zip(annotations, input.annotations))
        config_str = config["annotations"]["vep"]["params"]
        additional_str = "--vcf_info_field ANN --hgvsg"
        return f"{custom_str} {config_str} {additional_str}"
    else:
        return ""


rule get_vep_cache:
    output:
        directory("results/resources/vep/cache")
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    log:
        "logs/vep/cache.log"
    wrapper:
        "0.78.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("results/resources/vep/plugins")
    params:
        release=config["ref"]["release"]
    log:
        "logs/vep/plugins.log"
    wrapper:
        "0.78.0/bio/vep/plugins"


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]


def get_experiment_samples(experiment, case=True, control=True):
    case_samples, control_samples = [], []
    if case:
        case_samples = experiments[experiment]["samples"]["case"]
    if control:
        control_samples = experiments[experiment]["samples"]["control"]
    return case_samples + control_samples


def get_bam(wc):
    return samples.loc[wc.sample, "bam"]


def get_bai(wc):
    sample = samples.loc[wc.sample, "bam"]
    if sample.endswith(".cram"):
        return sample + ".crai"
    return sample + ".bai"


def get_vep_threads():
    n = len(samples)
    if n:
        return max(workflow.cores / n, 1)
    else:
        return 1


rule link_bai:
    input:
        "{x}.bai"
    output:
        "{x}.bam.bai"
    shell:
        "ln -sr {input} {output}"

rule link_crai:
    input:
        "{x}.crai"
    output:
        "{x}.cram.crai"
    shell:
        "ln -sr {input} {output}"


def get_vep_threads():
    n = len(samples)
    if n:
        return max(workflow.cores / n, 1)
    else:
        return 1


rule vcf_to_bcf_sample:
    input:
        vcf="results/strling/vcf/{experiment}/{sample}.vcf"
    output:
        bcf="results/strling/vcf/{experiment}/{sample}.bcf"
    shell:
        "bcftools sort -Ob {input.vcf} | bcftools norm -Ob -m-any > {output.bcf}"

rule vcf_to_bcf_ehdn_experiment:
    input:
        vcf="results/ehdn/final/{type}/{experiment}.vcf"
    output:
        bcf="results/ehdn/final/{type}/{experiment}.bcf"
    shell:
        "bcftools sort -Ob {input.vcf} | bcftools norm -Ob -m-any > {output.bcf}"
