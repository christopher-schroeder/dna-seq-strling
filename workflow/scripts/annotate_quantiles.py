import pysam
import numpy as np

with pysam.VariantFile(snakemake.input[0], "r") as f:
    header = f.header
    for i in [1,90,95,99]:
        header.add_meta('INFO', items=[('ID',f"Q{i}"), ('Number',1), ('Type','Float'), ('Description','{i:.2f} Quantile'.format(i=i/100))])
    with pysam.VariantFile(snakemake.output[0], "w", header=header) as o:
        for variant in f:
            controls = variant.info.get("CONTROLS", (0,))
            for i in [1,90,95,99]:
                variant.info[f"Q{i}"] = np.nanquantile(controls, i / 100)
            o.write(variant)