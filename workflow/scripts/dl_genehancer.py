import requests
import sys

with open(snakemake.output[0], "w") as o:
	fai=snakemake.input["fai"]
	hg=snakemake.params["hg"]

	for i, line in enumerate(open(fai, "r")):
		chrom, stop = line.split("\t")[0:2]
		stop = int(stop)
		if chrom.startswith("KI") or chrom.startswith("GL"):
			continue

		print(chrom, file=sys.stderr)

		data = {
			"hgsid": "0",
			"jsh_pageVertPos": "0",
			"clade": "mammal",
			"org": "Human",
			"db": f"{hg}",
			"hgta_group": "regulation",
			"hgta_track": "geneHancer",
			"hgta_table": "geneHancerRegElementsDoubleElite",
			"hgta_regionType": "range",
			"position": f"{chrom}:1-{stop:,}",
			"hgta_outputType": "primaryTable",
			"boolshad.sendToGalaxy": "0",
			"boolshad.sendToGreat": "0",
			"hgta_outFileName": "",
			"hgta_compressType": "none",
			"hgta_doTopSubmit": "get+output"
		}

		r = requests.post("https://genome.ucsc.edu/cgi-bin/hgTables", data=data)


		if i > 0:
			text = r.text.split("\n")
			text = "\n".join(text[1:])
		else:
			text=r.text

		print(text, file=o)