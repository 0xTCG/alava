#%%
import numpy as np
import statsmodels.api as sm
import glob
import sys
import pandas as pd
import math

pathway_path, magma_prefix = sys.argv[1:]
limit = 1.64

with open(pathway_path) as f:
    pathways = {(lp := l.strip().split())[0]: {"desc": lp[1]} for l in f}
    with open("_y.txt") as f:
        for l in f:
            l = l.strip().split()
            pathways["hsa" + l[0]]["genes"] = list(map(int, l[1:]))
print(f"Found {len(pathways)} pathways...")

zstat = {}
genes = None
for f in pathways.values():
    df = pd.read_fwf(f'{magma_prefix}.{f["desc"]}.glm.linear.magma.genes.out')
    if genes is None:
        genes = {int(r.GENE): (r.NPARAM, r.NSNPS) for _, r in df.iterrows()}
        print(f"Found {len(genes)} genes...")
    zstat[f["desc"]] = np.array(df.ZSTAT)

x0 = []
for g in genes:
    one_hot = [int(g in i["genes"]) for i in pathways.values()]
    NPARAM, NSNPS = genes[g]
    extra = [NPARAM, NPARAM / NSNPS, math.log(NPARAM), math.log(NPARAM / NSNPS)]
    x0.append(one_hot + extra)
x0 = np.array(x0)
x = sm.add_constant(x0)

#%%
print("Applying A-LAVA correction...")
for pathway, y in zstat.items():
    out_path = f'{magma_prefix}.{pathway}.glm.linear.alava'
    print(f"- {pathway} (file: {out_path})...")
    gls_model = sm.GLS(y, x)
    gls_results = gls_model.fit()
    results = []
    with open(out_path, 'w') as fo:
        for pi, p in enumerate(pathways):
            results.append((
                p,
                gls_results.params[pi + 1],
                gls_results.bse[pi + 1],
                gls_results.tvalues[pi + 1],
                gls_results.pvalues[pi + 1] / 2
            ))
            print(*results[-1][1:], file=fo)
    results.sort(key=lambda x: x[3], reverse=True)
    print("  Top calls:", ', '.join(f"{r[0]}: {float(r[3]):.3f}" for r in results if r[3] > limit))
