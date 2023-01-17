import yaml
import pandas as pd
from os.path import exists


data = pd.read_csv("samples.tsv", sep="\t").set_index("sample_name")
# for filename in data["bam"]:
#     if not exists(filename):
#         print(filename)

with open("experiments/expand.yaml", "r") as stream:
    try:
        experiment = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

for sample in experiment["samples"]["case"]:
    unit = data.loc[sample]
    filename = unit["bam"]
    if not exists(filename):
        print(filename)