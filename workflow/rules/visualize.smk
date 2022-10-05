rule visualize:
    input:
        bcf="results/{caller}/final/{type}/{experiment}.bcf"
    output:
        report(
            directory("results/{caller}/plots/{type}/{experiment}"),
            patterns=["{chrom}-{left}-{right}-{motif}-{gene}.pdf"],
            caption="../report/plots.rst",
            category="STR plots",
            subcategory="{experiment}",
        ),
    conda:
        "../envs/visualize.yaml"
    shell:
        "python workflow/scripts/visualize_single.py {input.bcf} {output}"