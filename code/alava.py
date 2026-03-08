#%%
import numpy as np
import statsmodels.api as sm
import glob
import sys
import pandas as pd
import math
import time
from math import log


LIMIT = 1.64


class DotDict(dict):
    def __getattr__(self, key):
        try: return self[key]
        except KeyError: raise AttributeError(key)
    def __setattr__(self, key, value): self[key] = value
    def __delattr__(self, key): del self[key]


class Timer:
    def __init__(self, name, header=True):
        self.name = name
        self.header = header

    def __enter__(self):
        self.start = time.perf_counter()
        if self.header: print(f"> {self.name}...")
        return self

    def __exit__(self, *args):
        self.end = time.perf_counter()
        print(f"> {self.name} took {self.end - self.start:.2f}s")


class Gene:
    def __init__(self, **kwargs):
        self.id = kwargs["GENE"]
        self.name = kwargs["ENTREZ"]
        self.chr = kwargs["CHR"]
        self.start = kwargs["START"]
        self.end = kwargs["STOP"]
        self.nsnp = kwargs["NSNPS"]
        self.nparam = kwargs["NPARAM"]
        self.zstat = kwargs["ZSTAT"]
        self.p = kwargs["P"]


class Pathway:
    def __init__(self, desc, genes):
        self.description = desc
        self.genes = genes
        self.binary = []  # binary indicator for genes
        self.zstat = []   # ZSTAT for each gene
        self.magma = {}
        self.alava = {}
        self.gene_info = {}

    def read_from_file(path):
        pathways = {}
        with open(path) as f:
            for l in f:
                l = l.strip().split()
                pathways[l[0]] = Pathway(l[1], list(map(int, l[2:])))
        return pathways

    def read_zstat(self, path, entrez_path=None):
        gene_ids = {}
        if entrez_path:
            with open(entrez_path) as f:
                for l in f:
                    l = l.strip().split()
                    if l[1].isdigit():
                        gene_ids[int(l[1])] = l[0].replace('"', '')
        df = pd.read_fwf(path)
        gene_info = {int(r.GENE): Gene(**{"ENTREZ": gene_ids.get(r.GENE, f"ID:{r.GENE}"), **r})
                     for _, r in df.iterrows()}
        self.binary = np.array([int(g in self.genes) for g in gene_info])
        self.zstat = np.array(df.ZSTAT)
        self.gene_info = gene_info
        return gene_info

    def calculate_magma(self, pathways, genes):
        self.magma = {}
        for pname, p in pathways.items():
            x0 = [[int(gn in p.genes), g.nparam, g.nparam / g.nsnp, log(g.nparam), log(g.nparam / g.nsnp)]
                  for gn, g in genes.items()]
            magma = sm.GLS(self.zstat, sm.add_constant(x0)).fit()
            self.magma[pname] = {
                "params": magma.params[1],
                "bse": magma.bse[1],
                "t": magma.tvalues[1],
                "p": magma.pvalues[1] / 2
            }

    def calculate_alava(self, pathways, genes):
        x0 = []
        for gn, g in genes.items():
            one_hot = [int(gn in p.genes) for p in pathways.values()]
            extra = [g.nparam, g.nparam / g.nsnp, log(g.nparam), log(g.nparam / g.nsnp)]
            x0.append(one_hot + extra)
        x = sm.add_constant(x0)

        alava = sm.GLS(self.zstat, x).fit()
        self.alava = {}
        for pi, pname in enumerate(pathways):
            self.alava[pname] = {
                "params": alava.params[pi + 1],
                "bse": alava.bse[pi + 1],
                "t": alava.tvalues[pi + 1],
                "p": alava.pvalues[pi + 1] / 2
            }

    def load_magma(self, path):
        with open(path) as f:
            for l in f:
                l = l.strip().split()
                self.magma[l[0]] = {
                    "params": float(l[1]), "bse": float(l[2]), "t": float(l[3]), "p": float(l[4]),
                }

    def load_alava(self, path):
        with open(path) as f:
            for l in f:
                l = l.strip().split()
                self.alava[l[0]] = {
                    "params": float(l[1]), "bse": float(l[2]), "t": float(l[3]), "p": float(l[4]),
                }

    @property
    def nice_name(self):
        return self.description[5:].replace("_", " ").title()


def print_results(results, fo=sys.stdout):
    results = sorted(
        [(rn, *r.values()) for rn, r in results.items()],
        key=lambda x: x[3], # t-score: results["t"]
        reverse=True
    )
    for r in results:
        print(*r, sep="\t", file=fo)


if __name__ == "__main__":
    pathway_path, magma_prefix = sys.argv[1:]

    with Timer("Reading pathway data"):
        pathways = Pathway.read_from_file(pathway_path)
        print(f"Found {len(pathways)} pathways...")

    gene_info = None
    with Timer("Reading gene and ZSTAT data"):
        for pname, p in pathways.items():
            gene_info = p.read_zstat(f'{magma_prefix}.{pname}.glm.linear.magma.genes.out')
        print(f"Found {len(gene_info)} genes...")

    with Timer("Applying MAGMA correction"):
        for pname, p in pathways.items():
            p.calculate_magma(pathways, gene_info)
            with open(f"{magma_prefix}.{pname}.glm.linear.magma", "w") as fo:
                print_results(p.magma, fo)

    with Timer("Applying ALAVA correction"):
        for pname, p in pathways.items():
            p.calculate_alava(pathways, gene_info)
            with open(f"{magma_prefix}.{pname}.glm.linear.alava", "w") as fo:
                print_results(p.alava, fo)


#%%

# x0 = []
# for g in genes:
#     one_hot = [int(g in i["genes"]) for i in pathways.values()]
#     NPARAM, NSNPS = genes[g]
#     extra = [NPARAM, NPARAM / NSNPS, math.log(NPARAM), math.log(NPARAM / NSNPS)]
#     x0.append(one_hot + extra)
# x0 = np.array(x0)
# x = sm.add_constant(x0)

#%%
# print("Applying MAGMA correction...")
# for pathway, y in zstat.items():
#     out_path = f'{magma_prefix}.{pathway}.glm.linear'
#     print(f"- Processing {pathway} (output: {out_path})...")

#     with open(out_path + ".magma", 'w') as fo:
#         for p in pathways:
#             magma_x0 = [
#                 [int(g in pathways[p]["genes"]),
#                 (NPARAM := genes[g][0]), NPARAM / (NSNPS := genes[g][1]),
#                 math.log(NPARAM), math.log(NPARAM / NSNPS)]
#                 for g in genes
#             ]
#             magma = sm.GLS(y, sm.add_constant(magma_x0)).fit()
#             result = {
#                 "name": p,
#                 "params": magma.params[1],
#                 "bse": magma.bse[1],
#                 "t": magma.tvalues[1],
#                 "p": magma.pvalues[1] / 2
#             }
#             print(*list(result.values()), sep="\t", file=fo)

# print("Applying A-LAVA correction...")
# for pathway, y in zstat.items():
#     out_path = f'{magma_prefix}.{pathway}.glm.linear'
#     print(f"- Processing {pathway} (output: {out_path})...")

#     alava = sm.GLS(y, x).fit()
#     results = []
#     with open(out_path + ".alava", 'w') as fo:
#         for pi, p in enumerate(pathways):
#             results.append({
#                 "name": p,
#                 "params": alava.params[pi + 1],
#                 "bse": alava.bse[pi + 1],
#                 "t": alava.tvalues[pi + 1],
#                 "p": alava.pvalues[pi + 1] / 2
#             })
#             print(*list(results[-1].values()), sep="\t", file=fo)
#     results.sort(key=lambda x: x["t"], reverse=True)
#     print("  Top A-LAVA calls:")
#     for r in results:
#         if r["t"] > limit:  # filter on t-value
#             print(f"    - {r['name']} ({pathways[r['name']]['desc']}): {float(r['t']):.3f}")
