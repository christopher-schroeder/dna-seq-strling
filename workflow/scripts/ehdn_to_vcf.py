from collections import namedtuple
import math
import argparse


#filename = "results/strling/outlier/1480-1.STRs.tsv"
#fileoutput = "/dev/stdout"


header = """##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=2,length=242193529>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=MT,length=16569>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##contig=<ID=KI270728.1,length=1872759>
##contig=<ID=KI270727.1,length=448248>
##contig=<ID=KI270442.1,length=392061>
##contig=<ID=KI270729.1,length=280839>
##contig=<ID=GL000225.1,length=211173>
##contig=<ID=KI270743.1,length=210658>
##contig=<ID=GL000008.2,length=209709>
##contig=<ID=GL000009.2,length=201709>
##contig=<ID=KI270747.1,length=198735>
##contig=<ID=KI270722.1,length=194050>
##contig=<ID=GL000194.1,length=191469>
##contig=<ID=KI270742.1,length=186739>
##contig=<ID=GL000205.2,length=185591>
##contig=<ID=GL000195.1,length=182896>
##contig=<ID=KI270736.1,length=181920>
##contig=<ID=KI270733.1,length=179772>
##contig=<ID=GL000224.1,length=179693>
##contig=<ID=GL000219.1,length=179198>
##contig=<ID=KI270719.1,length=176845>
##contig=<ID=GL000216.2,length=176608>
##contig=<ID=KI270712.1,length=176043>
##contig=<ID=KI270706.1,length=175055>
##contig=<ID=KI270725.1,length=172810>
##contig=<ID=KI270744.1,length=168472>
##contig=<ID=KI270734.1,length=165050>
##contig=<ID=GL000213.1,length=164239>
##contig=<ID=GL000220.1,length=161802>
##contig=<ID=KI270715.1,length=161471>
##contig=<ID=GL000218.1,length=161147>
##contig=<ID=KI270749.1,length=158759>
##contig=<ID=KI270741.1,length=157432>
##contig=<ID=GL000221.1,length=155397>
##contig=<ID=KI270716.1,length=153799>
##contig=<ID=KI270731.1,length=150754>
##contig=<ID=KI270751.1,length=150742>
##contig=<ID=KI270750.1,length=148850>
##contig=<ID=KI270519.1,length=138126>
##contig=<ID=GL000214.1,length=137718>
##contig=<ID=KI270708.1,length=127682>
##contig=<ID=KI270730.1,length=112551>
##contig=<ID=KI270438.1,length=112505>
##contig=<ID=KI270737.1,length=103838>
##contig=<ID=KI270721.1,length=100316>
##contig=<ID=KI270738.1,length=99375>
##contig=<ID=KI270748.1,length=93321>
##contig=<ID=KI270435.1,length=92983>
##contig=<ID=GL000208.1,length=92689>
##contig=<ID=KI270538.1,length=91309>
##contig=<ID=KI270756.1,length=79590>
##contig=<ID=KI270739.1,length=73985>
##contig=<ID=KI270757.1,length=71251>
##contig=<ID=KI270709.1,length=66860>
##contig=<ID=KI270746.1,length=66486>
##contig=<ID=KI270753.1,length=62944>
##contig=<ID=KI270589.1,length=44474>
##contig=<ID=KI270726.1,length=43739>
##contig=<ID=KI270735.1,length=42811>
##contig=<ID=KI270711.1,length=42210>
##contig=<ID=KI270745.1,length=41891>
##contig=<ID=KI270714.1,length=41717>
##contig=<ID=KI270732.1,length=41543>
##contig=<ID=KI270713.1,length=40745>
##contig=<ID=KI270754.1,length=40191>
##contig=<ID=KI270710.1,length=40176>
##contig=<ID=KI270717.1,length=40062>
##contig=<ID=KI270724.1,length=39555>
##contig=<ID=KI270720.1,length=39050>
##contig=<ID=KI270723.1,length=38115>
##contig=<ID=KI270718.1,length=38054>
##contig=<ID=KI270317.1,length=37690>
##contig=<ID=KI270740.1,length=37240>
##contig=<ID=KI270755.1,length=36723>
##contig=<ID=KI270707.1,length=32032>
##contig=<ID=KI270579.1,length=31033>
##contig=<ID=KI270752.1,length=27745>
##contig=<ID=KI270512.1,length=22689>
##contig=<ID=KI270322.1,length=21476>
##contig=<ID=GL000226.1,length=15008>
##contig=<ID=KI270311.1,length=12399>
##contig=<ID=KI270366.1,length=8320>
##contig=<ID=KI270511.1,length=8127>
##contig=<ID=KI270448.1,length=7992>
##contig=<ID=KI270521.1,length=7642>
##contig=<ID=KI270581.1,length=7046>
##contig=<ID=KI270582.1,length=6504>
##contig=<ID=KI270515.1,length=6361>
##contig=<ID=KI270588.1,length=6158>
##contig=<ID=KI270591.1,length=5796>
##contig=<ID=KI270522.1,length=5674>
##contig=<ID=KI270507.1,length=5353>
##contig=<ID=KI270590.1,length=4685>
##contig=<ID=KI270584.1,length=4513>
##contig=<ID=KI270320.1,length=4416>
##contig=<ID=KI270382.1,length=4215>
##contig=<ID=KI270468.1,length=4055>
##contig=<ID=KI270467.1,length=3920>
##contig=<ID=KI270362.1,length=3530>
##contig=<ID=KI270517.1,length=3253>
##contig=<ID=KI270593.1,length=3041>
##contig=<ID=KI270528.1,length=2983>
##contig=<ID=KI270587.1,length=2969>
##contig=<ID=KI270364.1,length=2855>
##contig=<ID=KI270371.1,length=2805>
##contig=<ID=KI270333.1,length=2699>
##contig=<ID=KI270374.1,length=2656>
##contig=<ID=KI270411.1,length=2646>
##contig=<ID=KI270414.1,length=2489>
##contig=<ID=KI270510.1,length=2415>
##contig=<ID=KI270390.1,length=2387>
##contig=<ID=KI270375.1,length=2378>
##contig=<ID=KI270420.1,length=2321>
##contig=<ID=KI270509.1,length=2318>
##contig=<ID=KI270315.1,length=2276>
##contig=<ID=KI270302.1,length=2274>
##contig=<ID=KI270518.1,length=2186>
##contig=<ID=KI270530.1,length=2168>
##contig=<ID=KI270304.1,length=2165>
##contig=<ID=KI270418.1,length=2145>
##contig=<ID=KI270424.1,length=2140>
##contig=<ID=KI270417.1,length=2043>
##contig=<ID=KI270508.1,length=1951>
##contig=<ID=KI270303.1,length=1942>
##contig=<ID=KI270381.1,length=1930>
##contig=<ID=KI270529.1,length=1899>
##contig=<ID=KI270425.1,length=1884>
##contig=<ID=KI270396.1,length=1880>
##contig=<ID=KI270363.1,length=1803>
##contig=<ID=KI270386.1,length=1788>
##contig=<ID=KI270465.1,length=1774>
##contig=<ID=KI270383.1,length=1750>
##contig=<ID=KI270384.1,length=1658>
##contig=<ID=KI270330.1,length=1652>
##contig=<ID=KI270372.1,length=1650>
##contig=<ID=KI270548.1,length=1599>
##contig=<ID=KI270580.1,length=1553>
##contig=<ID=KI270387.1,length=1537>
##contig=<ID=KI270391.1,length=1484>
##contig=<ID=KI270305.1,length=1472>
##contig=<ID=KI270373.1,length=1451>
##contig=<ID=KI270422.1,length=1445>
##contig=<ID=KI270316.1,length=1444>
##contig=<ID=KI270340.1,length=1428>
##contig=<ID=KI270338.1,length=1428>
##contig=<ID=KI270583.1,length=1400>
##contig=<ID=KI270334.1,length=1368>
##contig=<ID=KI270429.1,length=1361>
##contig=<ID=KI270393.1,length=1308>
##contig=<ID=KI270516.1,length=1300>
##contig=<ID=KI270389.1,length=1298>
##contig=<ID=KI270466.1,length=1233>
##contig=<ID=KI270388.1,length=1216>
##contig=<ID=KI270544.1,length=1202>
##contig=<ID=KI270310.1,length=1201>
##contig=<ID=KI270412.1,length=1179>
##contig=<ID=KI270395.1,length=1143>
##contig=<ID=KI270376.1,length=1136>
##contig=<ID=KI270337.1,length=1121>
##contig=<ID=KI270335.1,length=1048>
##contig=<ID=KI270378.1,length=1048>
##contig=<ID=KI270379.1,length=1045>
##contig=<ID=KI270329.1,length=1040>
##contig=<ID=KI270419.1,length=1029>
##contig=<ID=KI270336.1,length=1026>
##contig=<ID=KI270312.1,length=998>
##contig=<ID=KI270539.1,length=993>
##contig=<ID=KI270385.1,length=990>
##contig=<ID=KI270423.1,length=981>
##contig=<ID=KI270392.1,length=971>
##contig=<ID=KI270394.1,length=970>
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=END,Number=A,Type=Integer,Description="End position of structural variant (inclusive, 1-based).">
##INFO=<ID=SVTYPE,Number=A,Type=String,Description="Structural variant type">
##INFO=<ID=CONTROLS,Number=.,Type=Float,Description="The estimations of the controls.">
##FORMAT=<ID=A1,Number=A,Type=Float,Description="Estimation of the first allele.">
##FORMAT=<ID=A2,Number=A,Type=Float,Description="Estimation of the second allele.">
##FORMAT=<ID=QV,Number=A,Type=Float,Description="q-value.">"""

def convert(filename, fileoutput, cases):
    with open(filename, "r") as f, open(fileoutput, "w") as o:
        print(header, file=o)
        for i, line in enumerate(f):
            if i == 0:
                Line = namedtuple("Line", ", ". join(line.strip().split()))
                print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *cases, sep="\t", file=o)
                continue
            l = Line(*line.strip("\n").split())

            # if l.p_adj == "nan" or float(l.p_adj) > 0.1:
            #     continue

            counts = dict(map(lambda x: x.split(":"), l.counts.split(",")))

            # for key, value in counts.items():
            #     if 
            
            if float(l.bonf_pvalue) == 0.0:
                qual = 255
            elif float(l.bonf_pvalue) == 1.0:
                qual = 0.0
            else:
                qual = -10 * math.log10(float(l.bonf_pvalue))

            case_counts = "\t".join([f"{counts[c]}:.:." if c in counts else f".:.:." for c in cases])
            control_counts = ",".join([counts[c] for c in counts.keys() if c not in cases])
            
            if control_counts: 
                control_counts=f";CONTROLS={control_counts}"
            else:
                control_counts=""
    
            print(
                l.contig,
                l.start,
                ".",
                ".",
                l.motif,
                qual,
                "PASS",
                f"SVTYPE=STR;END={l.end}{control_counts}",
                "A1:A2:QV",
                case_counts,
                sep="\t", file=o)

filename = snakemake.input.tsv
fileoutput = snakemake.output.vcf
cases = snakemake.params.cases

# filename = "results/ehdn/casecontrol/aurora.casecontrol_locus.tsv"
# fileoutput = "/dev/stdout"
# cases = ["LNF-111_0", "LNF-117_0", "LNF-122_0", "LNF-122_1", "LNF-122_2", "LNF-131", "LNF-136_0", "LNF-136_1", "LNF-136_2", "LNF-35_0", "LNF-38_0", "LNF-44_0", "LNF-46_3", "LNF-55_0", "LNF-59_0", "LNF-61_0", "LNF-62_0", "LNF-65_0", "LNF-82_0", "LNF-82_1", "LNF-82_2", "LNF-99_0", "SPG-11_1", "SPG-142_0", "SPG-15_2", "SPG-212_0", "SPG-212_1", "SPG-212_2", "SPG-35_0", "SPG-36_0", "SPG-45_0", "SPG-59_0", "SPG-67_0", "SPG-68_0", "SPG-7_0", "SPG-7_1", "SPG-7_2", "SPG-77_0", "SPG-79_0", "SPG-94_0", "SPG-98_0"]
convert(filename, fileoutput, cases)