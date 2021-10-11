import glob

samples = pd.read_csv(config["samples"], sep="\t", comment='#', dtype={"sample_name": str, "bam": str, "group":str}).set_index("sample_name", drop=False).sort_index()
groups = samples["group"].unique()

groups_no_control = groups[groups != "control"]
samples_no_control = samples[samples["group"] != "control"]


wildcard_constraints:
    group="|".join(samples["group"].unique()),
    sample="|".join(samples["sample_name"]),


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
        "0.59.2/bio/vep/cache"


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]


def get_bam(wc):
    return samples.loc[wc.sample, "bam"]


def get_bai(wc):
    return samples.loc[wc.sample, "bam"] + ".bai"


def get_vep_threads():
    n = len(samples)
    if n:
        return max(workflow.cores / n, 1)
    else:
        return 1


rule get_vep_plugins:
    output:
        directory("results/resources/vep/plugins")
    params:
        release=config["ref"]["release"]
    log:
        "logs/vep/plugins.log"
    wrapper:
        "0.59.2/bio/vep/plugins"


rule link_bai:
    input:
        "{x}.bai"
    output:
        "{x}.bam.bai"
    shell:
        "ln -sr {input} {output}"