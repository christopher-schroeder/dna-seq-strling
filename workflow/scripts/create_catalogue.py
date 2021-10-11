import pandas as pd
import json

data = pd.read_csv(snakemake.input.regions, sep="\t")
ret = []

for i, row in data.iterrows():
    entry = dict()
    entry["LocusId"] = row["region_name"]
    entry["LocusStructure"] = f'({row["motif"]})*'
    entry["ReferenceRegion"] = f'{row["chrom"]}:{row["start"]}-{row["stop"]}'
    entry["VariantType"] = "Repeat"
    ret.append(entry)

with open(snakemake.output.catalogue, 'w') as f:
    json.dump(ret, f)