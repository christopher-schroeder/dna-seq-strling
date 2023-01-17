import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import argparse
import pysam
import math
import os

from collections import namedtuple
from matplotlib.lines import Line2D

import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("vcf", help="The STRs in vcf format.")
parser.add_argument("out")
args = parser.parse_args()


fs = 12  # fontsize

PlotData = namedtuple("PlotData", "chrom, left, right, motif, gene, consequence, values, samples")
Sample = namedtuple("Sample", "est1, est2, qv")
plotdata = []

os.makedirs(args.out, exist_ok=True)

with pysam.VariantFile(args.vcf) as f:
    sample_names = list(f.header.samples)

    for i, entry in enumerate(f):
        if len(entry.alts[0]) < 3:
            continue

        if not(entry.qual == None or entry.qual == "." or entry.qual > 0):
            continue

        if "CONTROLS" in entry.info:
            values = [z for x in entry.info["CONTROLS"] if not math.isnan(z:=float(x))]
        else:
            values = [0]

        # print(entry.info.keys())
        samples = []
        min_q_value = 1
        for j in range(len(entry.samples)):
            (est1,), (est2,), (qv,) = entry.samples[j].values()
            if qv is not None:
                min_q_value = min(min_q_value, qv)
            samples.append(Sample(est1, est2, qv))
        if min_q_value > 0.1 and entry.qual == ".":
            continue

        CSQ = entry.info["CSQ"][0]

        Allele, consequence, IMPACT, SYMBOL, Gene, Feature_type, Feature, BIOTYPE, EXON, INTRON, HGVSc, HGVSp, cDNA_position, CDS_position, Protein_position, Amino_acids, Codons, Existing_variation, DISTANCE, STRAND, FLAGS, SYMBOL_SOURCE, HGNC_ID = CSQ.split("|")
        consequence = ", ".join(map(lambda c: c.rsplit("_",1)[0], consequence.split("&")))
        plotdata.append(PlotData(entry.chrom, entry.pos, entry.stop, entry.alts[0], SYMBOL, consequence, values, tuple(samples)))

#plotdata = plotdata[0:18]

palette = sns.color_palette("colorblind", len(entry.samples))
plt.clf()

for i, d in enumerate(plotdata):
    fig, ax = plt.subplots()

    #sns.swarmplot(x=[1] * len(d.values), y=d.values, ax=ax, color="black", size=1.9)
    sns.violinplot(y=d.values, widths=0.3, showmeans=True, showextrema=True, showmedians=True, ax=ax, orient="v", linewidth=0, saturation=0.7, color="slategrey")

    draw_count = 0

    for j, s in enumerate(d.samples):
        drawn = False
        if s.est1 is not None:
            sns.lineplot(x=[-0.3, 0], y=[s.est1, s.est1], marker='o', markers=True, ax=ax, color=palette[j])
            plt.text(-0.35, s.est1, sample_names[j], horizontalalignment="right", color=palette[j])

        if s.est2 is not None:
            sns.lineplot(x=[0, 0.3], y=[s.est2, s.est2], marker='o', markers=True, ax=ax, color=palette[j], label=f"{float(s.qv):.2f}, {sample_names[j]}")
            plt.text(0.35, s.est2, sample_names[j], horizontalalignment="left", color=palette[j])

    ax.set_title(f'{d.chrom}:{d.left}-{d.right}, {d.motif}\n{d.consequence}, {d.gene}', fontsize=fs)
    ax.tick_params(axis='both', which='both', labelsize=10)
    ax.yaxis.set_tick_params(labelbottom=True)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax.set_xlim([-1, 1])


    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])

    plt.setp(ax.get_yticklabels(), visible=True)

    lines = [Line2D([0], [0], color=c, linewidth=1.5, linestyle='-') for c in palette]

    # legend1 = plt.legend(loc=2, title="q-values")

    #l4 = fig.legend(["1","2","3"], ncol=3 )
    # leg = fig.legend(lines, sample_names, loc="right", mode="expand")

    # plt.legend(lines, sample_names, loc='upper left', bbox_to_anchor=(1.0, 1.0), fancybox=False, shadow=False, ncol=1, title=f'{d.chrom}:{d.left}-{d.right}\n{d.motif}')

    # plt.subplots_adjust(hspace=0.4, left=0.05, right=0.95, top=1 - (0.4 / nrows), bottom=0.2 / nrows)
    # plt.margins(x=0, y=0)
    #leg.set_in_layout(True)
    # fig.tight_layout()

    # plt.gca().add_artist(legend1)
    fig.set_size_inches(10, 14)
    # plt.tight_layout()
    plt.margins(0.0, tight=True)
    plt.savefig(os.path.join(args.out, f"{d.chrom}-{d.left}-{d.right}-{d.motif}-{d.gene}.pdf"))
    plt.close()

    # print(i)
    # if i > 10:
    #     exit()