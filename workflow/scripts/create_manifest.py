import argparse
import os

# parser = argparse.ArgumentParser()
# parser.add_argument("--case", nargs="+")
# parser.add_argument("--control", nargs="+")
# args = parser.parse_args()
# case = args.case
# control = args.control

case = snakemake.input.case
control = snakemake.input.control

with open(snakemake.output.manifest, "w") as o:
    for filename in case:
        assert(filename.endswith(".str_profile.json"))
        sample_name = os.path.basename(filename[:-len(".str_profile.json")])
        print(sample_name, "case", filename, sep="\t", file=o)

    for filename in control:
        assert(filename.endswith(".str_profile.json"))
        sample_name = os.path.basename(filename[:-len(".str_profile.json")])
        print(sample_name, "control", filename, sep="\t", file=o)