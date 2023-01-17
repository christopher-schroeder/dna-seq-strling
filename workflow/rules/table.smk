def vembrane_expression(wc):
    samples = get_experiment_samples(wc.experiment, case=True, control=False)
    header = "index, chrom, pos, alt, region, gene, feature, distance, p-value, q1, q90, q95, q99"
    sample_header = expand("a1_{S}, a2_{S}, q_{S}", S=map(lambda s: s.replace("-", "_"), samples))
    header = ", ".join([header] + sample_header)

    expression = "INDEX, CHROM, POS, ALT, CSQ['Consequence'], CSQ['SYMBOL'], CSQ['Feature'], CSQ['DISTANCE'], 0, INFO['Q1'], INFO['Q90'], INFO['Q95'], INFO['Q99']" #10**(-QUAL/10)
    sample_expression = expand("FORMAT['A1']['{S}'], FORMAT['A2']['{S}'], FORMAT['QV']['{S}']", S=samples)
    expression = ", ".join([expression] + sample_expression)
    return header, expression


rule vembrane_table:
    input:
        bcf="results/{caller}/final/{type}/{experiment}.bcf"
    output:
        tsv=report(
            "results/{caller}/tables/{type}/{experiment}.tsv",
            caption="../report/table.rst",
            category="STR tables",
        )
    conda:
        "../envs/vembrane.yaml"
    params:
        expression=vembrane_expression
    log:
        "logs/vembrane-table/{caller}/{type}/{experiment}.log"
    shell:
        "bcftools norm -m-any {input} | vembrane table -k CSQ --header \"{params.expression[0]}\" \"{params.expression[1]}\" - > {output.tsv} 2> {log}"


rule tsv_to_excel:
    input:
        tsv="results/{x}.tsv"
    output:
        xlsx="results/{x}.xlsx"
    conda:
        "../envs/excel.yaml"
    log:
        "logs/tsv_to_xlsx/{x}.log"
    script:
        "../scripts/tsv_to_xlsx.py"
