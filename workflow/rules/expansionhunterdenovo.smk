rule ehdn_outlier:
    input:
        #case=lambda w: expand("results/profiles/{sample}.str_profile.json", sample=get_group_samples(w.group)),
        case=["results/ehdn/profiles/{sample}.str_profile.json"],
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=set(samples.sample_name) - set(get_group_samples(samples.loc[w.sample].group))),
        manifest="results/ehdn/manifests/{sample}.json",
        merged="results/ehdn/merged/{sample}.multisample_profile.json"
    output:
        "results/ehdn/outlier/{sample}.outlier_locus.tsv"
    log:
        "logs/ehdn/ehdn_casecontrol/{sample}.log"
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
        case=["results/ehdn/profiles/{sample}.str_profile.json"],
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=set(samples.sample_name) - set(get_group_samples(samples.loc[w.sample].group))),
        manifest="results/ehdn/manifests/{sample}.json",
        merged="results/ehdn/merged/{sample}.multisample_profile.json"
    output:
        "results/ehdn/casecontrol/{sample}.casecontrol_locus.tsv"
    log:
        "logs/strling/ehdn_casecontrol/{sample}.log"
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
        case=["results/ehdn/profiles/{sample}.str_profile.json"],
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=set(samples.sample_name) - set(get_group_samples(samples.loc[w.sample].group))),
        manifest="results/ehdn/manifests/{sample}.json",
        reference="resources/genome.fasta",
        fai="resources/genome.fasta.fai",
    output:
        merged="results/ehdn/merged/{sample}.multisample_profile.json"
    params:
        prefix="results/ehdn/merged/{sample}"
    conda:
        "../envs/expansionhunterdenovo.yaml"
    log:
        "logs/ehdn/merge/{sample}.log"
    shell:
        """ExpansionHunterDenovo merge \
        --manifest {input.manifest} \
        --reference {input.reference} \
        --output-prefix {params.prefix} 2> {log} 1>&2
        """


rule ehdn_create_manifest:
    input:
#        case=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=get_group_samples(w.group)),
        case=["results/ehdn/profiles/{sample}.str_profile.json"],
        control=lambda w: expand("results/ehdn/profiles/{sample}.str_profile.json", sample=set(samples.sample_name) - set(get_group_samples(samples.loc[w.sample].group))),
    output:
        manifest="results/ehdn/manifests/{sample}.json"
    conda:
        "../envs/expansionhunterdenovo.yaml"
    script:
        "../scripts/create_manifest.py"


rule ehdn_profile:
    input:
        bam=get_bam,
        bai=get_bai,
        reference="resources/genome.fasta"
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
