rule annotate_quantile_ehdn:
    input:
        "results/ehdn/vcf-experiment/outlier/{experiment}.bcf"
    output:
        "results/ehdn/vcf-experiment/outlier/{experiment}.quantile.bcf"
    script:
        "../scripts/annotate_quantiles_ehdn.py"


rule ehdn_outlier_to_vcf:
    input:
        tsv="results/ehdn/outlier/{experiment}.tsv"
    output:
        vcf="results/ehdn/vcf-experiment/{experiment}.vcf"
    params:
        cases=lambda w: get_experiment_samples(w.experiment, case=True, control=False)
    script:
        "../scripts/ehdn_outlier_to_vcf.py"


rule ehdn_outlier:
    input:
        case=   lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=True, control=False)),
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=False, control=True)),
        manifest="results/ehdn/manifests/{experiment}.json",
        merged="results/ehdn/merged/{experiment}.multisample_profile.json"
    output:
        "results/ehdn/outlier/{experiment}.tsv"
    log:
        "logs/ehdn/ehdn_outlier/{experiment}.log"
    shell:
        """
        python3 workflow/external/ehdn/scripts/outlier.py locus\
        --manifest {input.manifest} \
        --multisample-profile {input.merged} \
        --output {output} 2> {log} \
        """
    

rule ehdn_casecontrol:
    input:
        #case=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_group_samples(w.group)),
        case=   lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=True, control=False)),
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=False, control=True)),
        manifest="results/ehdn/manifests/{experiment}.json",
        merged="results/ehdn/merged/{experiment}.multisample_profile.json"
    output:
        "results/ehdn/casecontrol/{experiment}.tsv"
    log:
        "logs/strling/ehdn_casecontrol/{experiment}.log"
    shell:
        """
        python3 workflow/external/ehdn/scripts/casecontrol.py locus\
        --manifest {input.manifest} \
        --multisample-profile {input.merged} \
        --output {output} 2> {log} \
        """


rule ehdn_merge:
    input:
        #case=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_group_samples(w.group)),
        case=   lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=True, control=False)),
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=False, control=True)),
        manifest="results/ehdn/manifests/{experiment}.json",
        reference=genome,
        fai=genome + ".fai",
    output:
        merged="results/ehdn/merged/{experiment}.multisample_profile.json"
    params:
        prefix="results/ehdn/merged/{experiment}"
    conda:
        "../envs/expansionhunterdenovo.yaml"
    log:
        "logs/ehdn/merge/{experiment}.log"
    shell:
        """ExpansionHunterDenovo merge \
        --manifest {input.manifest} \
        --reference {input.reference} \
        --output-prefix {params.prefix} 2> {log} 1>&2
        """


rule ehdn_create_manifest:
    input:
#        case=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_group_samples(w.group)),
        case=   lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=True, control=False)),
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_experiment_samples(w.experiment, case=False, control=True)),
    output:
        manifest="results/ehdn/manifests/{experiment}.json"
    conda:
        "../envs/expansionhunterdenovo.yaml"
    script:
        "../scripts/create_manifest.py"


rule ehdn_profile:
    input:
        bam=get_bam,
        bai=get_bai,
        reference=genome
    output:
        "results/ehdn/profiles/{sample}.str_profile.json",
        "results/ehdn/profiles/{sample}.locus.tsv",
        "results/ehdn/profiles/{sample}.motif.tsv"
    params:
        prefix="results/ehdn/profiles/{sample}"
    conda:
        "../envs/expansionhunterdenovo.yaml"
    log:
        "log/ehdn/profile/{sample}.log"
    shell:
        """ExpansionHunterDenovo profile \
        --reads {input.bam} \
        --reference {input.reference} \
        --output-prefix {params.prefix} \
        --min-anchor-mapq 50 \
        --max-irr-mapq 40 2> {log}
        """
