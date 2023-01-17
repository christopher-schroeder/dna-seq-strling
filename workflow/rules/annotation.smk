rule annotate_vep:
    input:
        calls="results/{caller}/vcf-experiment/{type}/{experiment}.quantile.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls="results/{caller}/vcf-experiment/{type}/{experiment}.vep.bcf",
        stats="results/{caller}/vcf-experiment/{type}/{experiment}.vep.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra=""
    log:
        "logs/vep/{caller}/{type}/{experiment}.annotate.log"
    threads: get_vep_threads()
    wrapper:
        "0.78.0/bio/vep/annotate"


rule vembrane_filter:
    input:
        bcf="results/{caller}/vcf-experiment/{type}/{experiment}.vep.bcf",
    output:
        bcf="results/{caller}/vcf-experiment/{type}/{experiment}.vep.filtered.bcf",
    conda:
        "../envs/vembrane.yaml"
    params:
        expression="(QUAL and 10**(-QUAL/10) <= 0.1) or any(FORMAT['QV'][s] <= 0.1 for s in SAMPLES)"
    log:
        "logs/vembrane-table/{caller}/{type}/{experiment}.log"
    shell:
        "bcftools norm -m-any {input} | vembrane filter \"{params.expression}\" > {output}"


rule annotate_quantile:
    input:
        "results/{caller}/vcf-experiment/{type}/{experiment}.vep.filtered.bcf"
    output:
        "results/{caller}/vcf-experiment/{type}/{experiment}.quantile.bcf"
    script:
        "../scripts/annotate_quantiles.py"

rule just_copy:
    input:
        "results/{caller}/vcf-experiment/{type}/{experiment}.quantile.bcf",
    output:
        calls=report(
            "results/{caller}/final/{type}/{experiment}.bcf",
            caption="../report/str.rst",
            category="Test",
        ),
    shell:
        "cat {input} > {output}"


# rule annotate_genehancer:
#     input:
#         config="workflow/resources/genehancer.yaml",
#         bcf="results/{caller}/vcf-experiment/{type}/{experiment}.vep.bcf",
#         genehancer="results/resources/genehancer.tsv"
#     output:
#         calls=report(
#             "results/{caller}/final/{type}/{experiment}.bcf",
#             caption="../report/str.rst",
#             category="Test",
#         ),
#     log:
#         "logs/annotate_genehancer/{caller}/{type}/{experiment}.log"
#     conda:
#         "../envs/vembrane.yaml"
#     shell:
#         "vembrane annotate {input.config} {input.bcf} > {output} 2> {log}"