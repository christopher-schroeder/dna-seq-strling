rule annotate_vep:
    input:
        calls="results/strling/vcf/{experiment}/{experiment}.all.q.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls=pipe("results/strling/vcf/{experiment}/{experiment}.all.vep.bcf"),
        stats="results/strling/vcf/{experiment}/{experiment}.all.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra=""
    log:
        "logs/vep/{experiment}.annotate.log"
    threads: get_vep_threads()
    wrapper:
        "0.78.0/bio/vep/annotate"


rule annotate_genehancer:
    input:
        config="workflow/resources/genehancer.yaml",
        bcf="results/strling/vcf/{experiment}/{experiment}.all.vep.bcf",
        genehancer="results/resources/genehancer.tsv"
    output:
        calls=report(
            "results/strling/vcf/{experiment}/{experiment}.all.annotated.bcf",
            caption="../report/str.rst",
            category="Test",
        ),
    log:
        "logs/annotate_genehancer/{experiment}.log"
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane annotate {input.config} {input.bcf} > {output} 2> {log}"