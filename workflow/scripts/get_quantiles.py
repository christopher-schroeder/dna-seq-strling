import pandas as pd
import numpy as np
import argparse
import math

from collections import defaultdict

# parser = argparse.ArgumentParser()
# parser.add_argument("genotypes", nargs="*")
# args = parser.parse_args()

# filenames = args.genotypes
fileoutput = "/dev/stdout"

#filenames = [
#    "/vol/nano/depienne/strling/pipeline/results/strling/call/Homo_219-genotype.txt",
#    "/vol/nano/depienne/strling/pipeline/results/strling/call/Homo_247-genotype.txt",
#]

filenames = snakemake.input.genotypes
fileoutput = snakemake.output.quantiles

allele = defaultdict(lambda: defaultdict(list))
with open(fileoutput, "w") as o:
    for f in filenames:
        data = pd.read_csv(f, sep="\t")
        allele_est = data[["#chrom", "left", "right", "repeatunit", "allele1_est", "allele2_est"]]

        for chrom, left, right, repeatunit, allele1_est, allele2_est in allele_est.values:
            allele[(chrom, left, right, repeatunit)][0].append(allele1_est)
            allele[(chrom, left, right, repeatunit)][1].append(allele2_est)

    for key, value in allele.items():
        values_1 = [v for v in value[0] if not math.isnan(v)]
        values_2 = [v for v in value[1] if not math.isnan(v)]
        
        # if not values_1:
        #     values_1 = [float("nan")]
        # if not values_2:
        #     values_2 = [float("nan")]

        values = values_1 + values_2

        if len(values) == 0:
            values = [float("nan")]

        #quantile_alelle1 = np.quantile(values_1, [0.01, 0.99])
        #quantile_alelle2 = np.quantile(values_2, [0.01, 0.99])
        #print(*key, *quantile_alelle1, *quantile_alelle2, sep="\t")

        quantile_alelle = np.quantile(values, [0.01 ,0.5, 0.9, 0.95, 0.99])

        print(*key, *quantile_alelle, sep="\t", file=o)