import pandas as pd
from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import math
import sys

from collections import defaultdict, OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("genotypes", nargs="*")
args = parser.parse_args()

filenames = args.genotypes
fileoutput = "/dev/stdout"

#filenames = [
#    "/vol/nano/depienne/strling/pipeline/results/strling/call/Homo_219-genotype.txt",
#    "/vol/nano/depienne/strling/pipeline/results/strling/call/Homo_247-genotype.txt",
#]

# filenames = snakemake.input.genotypes
# fileoutput = snakemake.output.quantiles
    
d_values = dict()

allele = defaultdict(list)
with open(fileoutput, "w") as o:
    for f in filenames:
        data = pd.read_csv(f, sep="\t")
        allele_est = data[["#chrom", "left", "right", "repeatunit", "allele1_est", "allele2_est"]]

        for chrom, left, right, repeatunit, allele1_est, allele2_est in allele_est.values:
            allele[(chrom, left, right, repeatunit)].append(allele1_est)
            allele[(chrom, left, right, repeatunit)].append(allele2_est)


last_line = None

for line in sys.stdin:
    if line.startswith("#"):
        if line.startswith("##FORMAT") and not last_line.startswith("##FORMAT"):
            print('##INFO=<ID=CONTROLS,Number=.,Type=Float,Description="The estimations of the controls.">')
        print(line, end="")
        last_line = line
        continue
    
    split = line.split("\t")
    chrom, left = split[0:2]
    motif = split[4]
    info = OrderedDict(x.split("=") for x in split[7].split(";"))
    right = info["END"]

    key = (chrom, int(left), int(right), motif)
    
    values = allele[key]

    if len(values) == 0:
        values = [float("nan")]


    info["CONTROLS"] = ",".join(map(str, values))
    split[7] = ";".join(f"{k}={v}" for k,v in info.items())

    print(*split, sep="\t", end="")