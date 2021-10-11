import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import argparse
import pysam
import math
from collections import namedtuple
from matplotlib.lines import Line2D

import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("out")
args = parser.parse_args()


fs = 12  # fontsize

PlotData = namedtuple("PlotData", "chrom, left, right, motif, values, samples")
Sample = namedtuple("Sample", "est1, est2, qv")
plotdata = []

with pysam.VariantFile(args.vcf) as f:
    sample_names = list(f.header.samples)

    for i, entry in enumerate(f):
        if len(entry.alts[0]) < 3:
            continue

        values = [z for x in entry.info["CONTROLS"] if not math.isnan(z:=float(x))]

        samples = []
        min_q_value = 1
        for j in range(len(entry.samples)):
            (est1,), (est2,), (qv,) = entry.samples[j].values()
            if qv is not None:
                min_q_value = min(min_q_value, qv)
            samples.append(Sample(est1, est2, qv))

        if min_q_value > 0.1:
            continue

        plotdata.append(PlotData(entry.chrom, entry.pos, entry.stop, entry.alts[0], values, tuple(samples)))

#plotdata = plotdata[0:18]

palette = sns.color_palette("colorblind", len(entry.samples))

ncols = 6
nrows = math.ceil(len(plotdata) / ncols)

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 3, nrows * 3))

for i, d in enumerate(plotdata):
    x, y = i // ncols, i % ncols
    if nrows == 1:
        ax = axes[y]
    else:
        ax = axes[x, y]

    #sns.swarmplot(x=[1] * len(d.values), y=d.values, ax=ax, color="black", size=1.9)
    sns.violinplot(y=d.values, widths=0.3, showmeans=True, showextrema=True, showmedians=True, ax=ax, orient="v", linewidth=0, saturation=0.7, color="slategrey")

    for j, s in enumerate(d.samples):
        if s.est1 is not None:
            sns.lineplot(x=[-1, 0], y=[s.est1, s.est1], marker='o', markers=True, ax=ax, color=palette[j])
        if s.est2 is not None:
            sns.lineplot(x=[0, 1], y=[s.est2, s.est2], marker='o', markers=True, ax=ax, color=palette[j], label="{:.2f}".format(float(s.qv)))

    ax.set_title(f'{d.chrom}:{d.left}-{d.right}\n{d.motif}', fontsize=fs)
    ax.tick_params(axis='both', which='both', labelsize=10)
    ax.yaxis.set_tick_params(labelbottom=True)
    ax.legend(loc="upper left")
    ax.set_xlim([-1, 1])

    plt.setp(ax.get_yticklabels(), visible=True)

# i+=1
# while i < ncols * nrows:
#     x, y = i // ncols, i % ncols
#     ax.set_visible(False)
#     i += 1

lines = [Line2D([0], [0], color=c, linewidth=1.5, linestyle='-') for c in palette]

#l4 = fig.legend(["1","2","3"], ncol=3 )
leg = fig.legend(lines, sample_names, loc="upper center", mode="expand", ncol=min(6, len(sample_names)))
plt.subplots_adjust(hspace=0.4, left=0.05, right=0.95, top=1 - (0.4 / nrows), bottom=0.2 / nrows)
#plt.margins(x=0, y=0)
#leg.set_in_layout(True)
#fig.tight_layout()
plt.savefig(args.out)

# 0.6