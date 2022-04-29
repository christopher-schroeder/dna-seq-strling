import glob
import pandas as pd

samples = pd.read_csv(config["samples"], sep="\t", comment='#', dtype={"sample_name": str, "bam": str, "group":str}).set_index("sample_name", drop=False).sort_index()
groups = samples["group"].unique()

groups_no_control = groups[groups != "control"]
samples_no_control = samples[samples["group"] != "control"]
experiments = config["experiments"]

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
    groups = []
    if case:
        groups += config["experiments"][experiment]["case"]
    if control:
        groups += config["experiments"][experiment]["control"]
    return samples.loc[samples["group"].isin(groups)]["sample_name"]


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


def get_vep_threads():
    n = len(samples)
    if n:
        return max(workflow.cores / n, 1)
    else:
        return 1