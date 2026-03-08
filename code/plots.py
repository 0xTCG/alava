#%%
import os, glob, re
os.environ['R_HOME'] = '/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/lib64/R'

import numpy as np
import pandas as pd
import alava
import tabulate
from plotnine import *
from scipy.stats import pearsonr
import collections
import sys
import math
from rpy2.robjects import r, NULL
from rpy2.robjects.packages import importr
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri


# Parameters: T-value and p-value thresholds
T, P = 1.64, 0.05
RESULTS, PREFIX, INPUT, PLOTS = sys.argv[1:]
# RESULTS = "results/out"
# PREFIX = f"{RESULTS}/eur/magma/hg19-qc"
# INPUT = "data/simulations.txt"
# PLOTS = f"{RESULTS}/plots"

# PREFIX, PLOTS = "out-old/all_chr.rsq0.5_EUR_QC_annotated", "out-old/plots"

# Read gene data
pathways = alava.Pathway.read_from_file(f"{RESULTS}/kegg-pathways.txt")
for pname, p in pathways.items():
    gene_info = p.read_zstat(f"{PREFIX}.{pname}.glm.linear.magma.genes.out",
                             f"{RESULTS}/entrez-ids.txt")

#%% Figure 1: 📊 use as-is
### Figure 2: 📊 use as-is
### Figure 3, 4: 🚧 see sim.py

#%% Figure 5, 6:
# Interesting pathways
interesting = [ "hsa00982", "hsa00030", "hsa00920", ]
counts = []
for pathway in interesting:
    r = pathways[pathway]
    r.load_alava(f"{PREFIX}.{pathway}.glm.linear.alava")
    r.load_magma(f"{PREFIX}.{pathway}.glm.linear.magma")

    alava = sorted([(p, v) for p, v in r.alava.items() if v["t"] >= T], key=lambda x: x[1]["p"], reverse=False)
    magma = sorted([(p, v) for p, v in r.magma.items() if v["t"] >= T], key=lambda x: x[1]["p"], reverse=False)
    print(f"{r.description} limits:", len(alava), len(magma))
    table = []
    for (m, md), (a, ad) in zip(alava, magma):
        if not ad["p"] <= P and not md["p"] <= P:
            break
        if ad["p"] <= P:
            counts.append(a)
        table.append(
            [f"{m} ({pathways[m].nice_name})", md["params"], md["p"],
             f"{a} ({pathways[a].nice_name})", ad["params"], ad["p"]]
        )
    with open(f"{PLOTS}/f5_{pathway}.txt", "w") as fo:
        print(tabulate.tabulate(table, headers=["Pathway", "ß", "p"] * 2), file=fo)

    alava = pd.DataFrame(r.alava).transpose()
    magma = pd.DataFrame(r.magma).transpose()
    df = alava.reset_index().merge(magma.reset_index(), on="index", suffixes=("_alava", "_magma"))
    df['logp_alava'] = -np.log10(df['p_alava'])
    df['logp_magma'] = -np.log10(df['p_magma'])

    pp = (
        ggplot(df, aes(x='logp_alava', y='logp_magma'))
        + geom_point(alpha=0.5, color="red")
        + geom_abline(intercept=0, slope=1, linetype='dashed')
        + annotate(
            "text", x=0.8, y=3.5,
            label=f"Pearson r: {pearsonr(df['logp_alava'], df['logp_magma'])[0]:.3f}"
        )
        + xlim(0, 3.5) + ylim(0, 3.5)
        + labs(title=r.nice_name,
               subtitle=f"-log(p) correlation", x="ALAVA", y="MAGMA")
    )
    pb = (
        ggplot(df, aes(x='params_alava', y='params_magma'))
        + geom_point(alpha=0.5, color="green")
        + geom_abline(intercept=0, slope=1, linetype='dashed')
        + annotate(
            "text", x=-0.55, y=1,
            label=f"Pearson r: {pearsonr(df['params_alava'], df['params_magma'])[0]:.3f}"
        )
        + xlim(-1, 1) + ylim(-1, 1)
        + labs(subtitle=f"ß-value correlation", x="ALAVA", y="")
    )
    (
        (pp | pb)
        & theme_light()
        & theme(plot_title_position="plot", aspect_ratio=1)
    ).save(filename=f"{PLOTS}/f6_{pathway}.png")

#%% Appendix: Figure 7
p = []
for pathway in pathways:
    r = pathways[pathway]
    for g in sorted([g for g in r.gene_info.values() if g.p <= 0.0003], key=lambda x: x.p):
        p.append(g.name)
with open(f"{PLOTS}/f7.txt", "w") as fo:
    for n, c in collections.Counter(p).most_common():
        print(c, n, sep="\t", file=fo)

### Appendix: Figure 8 ❌
### Appendix: Figure 9, 10 🚧 see Figure 6

#%% Appendix: Figure 11, 12, 14:
### Appendix: Figure 13: 🤔 (use as-is?)
genes_df = None
for fname in glob.glob(f"{PREFIX}.*.genes.out"):
    d = pd.read_fwf(fname)
    d = d[d.P < P]
    d["KEGG"] = re.match(r".+\.(hsa\d+)\..+", fname)[1]
    d["ENTREZ"] = d.GENE.apply(lambda x: gene_info[x].name)
    if genes_df is None: genes_df = d
    else: genes_df = pd.concat([genes_df, d])
data_table = importr("data.table")
with localconverter(default_converter + pandas2ri.converter):
    r_df = pandas2ri.py2rpy(genes_df)
r_dt = data_table.as_data_table(r_df)

# Import R libraries
qqman = importr("qqman")
r(f'png("{PLOTS}/f11.png", width=1600, height=900)')
qqman.manhattan(
    r_dt,
    chr="CHR",
    bp="START",
    p="P",
    snp="ENTREZ",
    col=r.c("darkcyan","darkblue","skyblue","lightblue","lightgray"),
    chrlabs=NULL,
    suggestiveline=float(-math.log10(2.9e-05)),
    genomewideline=float(-math.log10(2.8e-06)),
    highlight=NULL,
    logp=True,
    annotatePval=0.000028,
    annotateTop=False
)
r("dev.off()")

snps_df = None
for fname in glob.glob(f"{PREFIX}.*.linear"):
    d = pd.read_table(fname)
    d = d[d.P < 0.01]
    d["KEGG"] = re.match(r".+\.(hsa\d+)\..+", fname)[1]
    if snps_df is None: snps_df = d
    else: snps_df = pd.concat([snps_df, d])
    print(fname, len(snps_df))
with localconverter(default_converter + pandas2ri.converter):
    r_df = pandas2ri.py2rpy(snps_df)
r_dt = data_table.as_data_table(r_df)

r(f'png("{PLOTS}/f12.png", width=1200, height=1200)')
qqman.qq(r_dt.rx2("P"))
r("dev.off()")

GWASinspector = importr("GWASinspector")
GWASinspector.manhattan_plot(
    r_dt,
    chr="#CHROM",
    pvalue="P",
    position="POS",
    fileName=f"{PLOTS}/f14.png",
    plot_title="SNP-level Manhattan Plot",
    plot_subtitle="",
    p_threshold=0.9,
    sig_threshold_log=float(-math.log10(5 * 10**-8)),
    beta="ID",
    std_error=NULL,
    check_columns=True
)

#%% Appendix: Figure 15:
table = []
for pathway in pathways:
    r = pathways[pathway]
    for gn, g in r.gene_info.items():
        if g.p <= 2.9e-5:
            table.append((g.p, g.name, g.chr, r.nice_name))
with open(f"{PLOTS}/f15.txt", "w") as fo:
    print(tabulate.tabulate(
            sorted(table, key=lambda x: x[0]), headers=["p-value", "Gene", "Chromosome", "Trait"]
          ), file=fo)

#%% Appendix: Figure 16 & 17:
table = []
for p, cnt in collections.Counter(counts).most_common():
    mean_p = np.average([r.alava[p]["p"] for r in pathways.values() if p in r.alava])
    mean_b = np.average([r.alava[p]["params"] for r in pathways.values() if p in r.alava])
    table.append((cnt, p, pathways[p].nice_name, mean_p, mean_b))
with open(f"{PLOTS}/f16.txt", "w") as fo:
    print(tabulate.tabulate(table, headers=["Count", "KEGG ID", "Trait", "p", "ß"]),
          file=fo)

### Appendix: Figure 18, 19, 20, 21 ❌
### Appendix: Figure 22, 23, 24, 25, 26, 27  🚧 see R code in plots.ipynb
### Appendix: Figure 28, 29, 30  🚧 see Figure 15
### Appendix: Figure 31 📊 use as-is


#%% Extra: Mansoureh's original KEGG list
if False:
    pathways = {}
    with (
        open("/scratch/inumanag/alava/final-instr/sim/Part1/SIMULATION_DATA_FILES/PREDEFINED_INFORMATION/65_metabolic_KEGG_Ids_2") as f1,
        open("/scratch/inumanag/alava/final-instr/sim/Part1/SIMULATION_DATA_FILES/PREDEFINED_INFORMATION/65_hsa_KEGG_pathways_genes_reference") as f2
    ):
        for l in f1:
            l = l.strip().split()
            pathways[l[0]] = alava.Pathway(l[1], [])
        for l in f2:
            l = l.strip().split()
            pathways["hsa" + l[0]].genes = list(map(int, l[1:]))
    with open("out-old/kegg-pathways.txt", "w") as fo:
        for p in sorted(pathways):
            print(p, pathways[p].description, ' '.join(str(x) for x in pathways[p].genes), file=fo)

#%% Extra: Rename Mansoureh -> current
if False:
    import re, shutil
    for g in glob.glob("out-old/*.KEGG_*"):
        if (m := re.match(r".+\.(KEGG_[A-Z0-9_]+)\.", g)):
            i = [i for i in pathways if pathways[i].description == m[1]][0]
            nf = g.replace(m[1], i)
            if not os.path.exists(nf):
                shutil.move(g, nf)
