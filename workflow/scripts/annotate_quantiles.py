import pandas as pd
from collections import defaultdict

# quantiles_filename = "results/strling/quantiles.txt"
# filename = "results/strling/outlier/ataxia_S012.STRs.tsv"

filename = snakemake.input.outliers
quantiles_filename = snakemake.input.quantiles
fileoutput = snakemake.output[0]

quantiles_data = pd.read_csv(quantiles_filename, sep="\t", header=None).iterrows()
quantiles = {(chrom, left, right, motif):quantiles for index, (chrom, left, right, motif, *quantiles) in quantiles_data}
quantiles = defaultdict(lambda: [float("nan"), float("nan"), float("nan"), float("nan"), float("nan")], quantiles)

with open(fileoutput, "w") as o:
    data = pd.read_csv(filename, sep="\t", header=0)
    print(*data.keys().values, *["q1", "q5", "q90", "q95", "q99"], sep="\t", file=o)
    for index, row in data.iterrows():
        key = tuple(row[["chrom", "left", "right", "repeatunit"]].values)
        print(*row, *quantiles[key], sep="\t", file=o)