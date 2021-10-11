rule annotate_variants:
    input:
        calls="results/strling/vcf/{group}/{group}.all.q.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls=report(
            "results/strling/vcf/{group}/{group}.all.annotated.bcf",
            caption="../report/str.rst",
            category="Test",
            ),
        stats="results/strling/vcf/{group}/{group}.all.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra=""
    log:
        "logs/vep/{group}.annotate.log"
    threads: get_vep_threads()
    wrapper:
        "0.78.0/bio/vep/annotate"