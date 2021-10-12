def vembrane_expression(wc):
    samples = get_group_samples(wc.group)
    header = "index, chrom, pos, alt, region, gene, feature, distance"
    sample_header = expand("a1_{S}, a2_{S}, q_{S}", S=map(lambda s: s.replace("-", "_"), samples))
    header = ", ".join([header] + sample_header)

    expression = "INDEX, CHROM, POS, ALT[0], CSQ['Consequence'], CSQ['SYMBOL'], CSQ['Feature'], CSQ['DISTANCE']"
    sample_expression = expand("FORMAT['A1']['{S}'], FORMAT['A2']['{S}'], FORMAT['QV']['{S}']", S=samples)
    expression = ", ".join([expression] + sample_expression)

    return header, expression


rule vembrane_table:
    input:
        bcf="results/strling/vcf/{group}/{group}.all.annotated.filtered.bcf",
    output:
        tsv=report(
            "results/strling/tables/{group}/{group}.tsv",
            caption="../report/table.rst",
            category="{group}",
        )
    conda:
        "../envs/vembrane.yaml"
    params:
        expression=vembrane_expression
    log:
        "logs/vembrane-table/{group}.log"
    shell:
        "bcftools norm -m-any {input} | vembrane table -k CSQ --header \"{params.expression[0]}\" \"{params.expression[1]}\" - > {output.tsv} 2> {log}"


rule vembrane_filter:
    input:
        bcf="results/strling/vcf/{group}/{group}.all.annotated.bcf",
    output:
        bcf="results/strling/vcf/{group}/{group}.all.annotated.filtered.bcf",
    conda:
        "../envs/vembrane.yaml"
    params:
        expression="any(FORMAT['QV'][s] <= 0.1 for s in SAMPLES)"
    log:
        "logs/vembrane-table/{group}.log"
    shell:
        "bcftools norm -m-any {input} | vembrane filter \"{params.expression}\" > {output}"


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
