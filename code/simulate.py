#%%
import sys
import random
import numpy as np
import alava
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as mpatches
import tabulate


RESULTS, PREFIX, INPUT, PLOTS = sys.argv[1:]
# RESULTS = "results/out"
# PREFIX = f"{RESULTS}/eur/magma/hg19-qc"
# INPUT = "data/simulations.txt"
# PLOTS = f"{RESULTS}/plots"
random.seed(5)


#%%
pathways = alava.Pathway.read_from_file(f"{RESULTS}/kegg-pathways.txt")
for pname, p in pathways.items():
    gene_info = p.read_zstat(f"{PREFIX}.{pname}.glm.linear.magma.genes.out",
                             f"{RESULTS}/entrez-ids.txt")

simulations = {}
with open(INPUT) as f:
    for l in f:
        name, overlap_sets = l.split(": ")
        name = f"hsa{name}"
        simulations[name] = alava.DotDict()
        simulations[name].overlaps = [f"hsa{i}" for i in overlap_sets.strip().split()]
        assert name in pathways
        assert all(i in pathways for i in simulations[name].overlaps)

def generate_random_simulation(overlaps, fn):
    frequency_sig = sum(pathways[p].binary for p in overlaps)
    for v in frequency_sig:
        rnd = fn(v)
        yield rnd

fns = {
    0: (lambda v: random.uniform(0, 1.5) if v == 0 else random.uniform(2.5, 3)),
    # Other generator functions that can be tried:
    # 1: (lambda v: random.uniform(0, 1.5) if v == 0 else random.uniform(3, 5.5)),
    # 2: (lambda v: 0 if v == 0 else 3),
    # 3: (lambda v: 0 if v == 0 else max(3 + v, 6)),
}

def generate_simulation(data):
    overlaps, pathways, pathway = data
    experiments = {}
    for fn in fns:
        pathways[pathway].zstat = np.array(np.fromiter(
            generate_random_simulation(overlaps, fns[fn]), dtype=np.float64
        ))
        pathways[pathway].calculate_magma(pathways, gene_info)
        pathways[pathway].calculate_alava(pathways, gene_info)
        experiments[fn] = alava.DotDict({
            "zstat": pathways[pathway].zstat,
            "magma": pathways[pathway].magma,
            "alava": pathways[pathway].alava,
        })
    return (pathway, experiments)

for p, s in simulations.items():
    with alava.Timer(f"Simulation {p}", False):
        s.experiments = generate_simulation((s.overlaps, pathways, p))
print(f"Processed {len(simulations)} simulations")

#%%
confusion = {"TN": 0, "FP": 1, "TP": 2, "FN": 3}
index = {n: i for i, n in enumerate(simulations)}

def plots(fn):
    mats = {
        "true": [[confusion["TP"] if s in ss.overlaps else confusion["TN"] for s in simulations]
                for ss in simulations.values()],
        "magma": [[confusion["TP"] if s in ss.overlaps else confusion["TN"] for s in simulations]
                for ss in simulations.values()],
        "alava": [[confusion["TP"] if s in ss.overlaps else confusion["TN"] for s in simulations]
                for ss in simulations.values()]
    }

    total_genes, true_sig, alava_sig, magma_sig = [], [], [], []
    results = {}
    for pname, simulation in simulations.items():
        results[pname] = {"true": sorted(simulations[pname].overlaps)}
        data = simulation.experiments[1][fn]
        for method in ("alava", "magma"):
            results[pname][method] = []
            for p in pathways:
                b, pv = round(data[method][p]["params"], 2), round(data[method][p]["p"], 2)
                if (found := b >= 0 and pv <= 0.05):
                    results[pname][method].append(p)

                if mats[method][index[pname]][index[p]] == confusion["TN"] and found:
                    mats[method][index[pname]][index[p]] = confusion["FP"]
                elif mats[method][index[pname]][index[p]] == confusion["TP"] and not found:
                    mats[method][index[pname]][index[p]] = confusion["FN"]
        total_genes.append(sum(1 for i in data.zstat if i >= 1.5))  # Number of significant genes
        true_sig.append(len(results[pname]["true"]))
        alava_sig.append(len(results[pname]["alava"]))
        magma_sig.append(len(results[pname]["magma"]))

    plt.scatter(total_genes, true_sig, label="Expected")
    plt.scatter(total_genes, alava_sig, label="A-LAVA")
    plt.scatter(total_genes, magma_sig, label = "MAGMA")
    plt.xlabel("Number of Significant Genes")
    plt.ylabel("Number of Calculated Significant Pathways")
    plt.legend(bbox_to_anchor=(1, 1))
    plt.savefig(f"{PLOTS}/fig3.png")

    cols = ["#ffffff", "#f02e35", "#8fe56f", "#981534"]
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7, 3), layout='constrained', squeeze=False)
    images = []
    for [ax, data] in zip(axs.flat, [mats["alava"], mats["magma"]]):
        images.append(ax.pcolormesh(data, cmap=colors.LinearSegmentedColormap.from_list("", cols)))
    axs[0, 0].set_title("A-LAVA")
    axs[0, 1].set_title("MAGMA")
    fig.legend(handles=[mpatches.Patch(color=cols[i], label=label) for i, label in enumerate(confusion)],
            loc='outside right')
    plt.savefig(f"{PLOTS}/fig4.png")

    table = {}
    for method in ("true", "alava", "magma"):
        D = {}
        for i, label in enumerate(confusion):
            D[label] = sum(1 for r in mats[method] for c in r if c == i)
            table.setdefault(method, []).append(D[label])
        table[method].append(round(P := float(D["TP"]) / (D["TP"] + D["FP"]), 3))
        table[method].append(round(R := float(D["TP"]) / (D["TP"] + D["FN"]), 3))
        table[method].append(round(2 * P * R / (P + R), 3))

    with open(f"{PLOTS}/fig4.txt", "w") as fo:
        print(
            tabulate.tabulate(
                table, headers="keys",
                showindex=list(confusion.keys()) + ["Precs.", "Recall", "F1"],
                numalign="right",
                tablefmt="rst",
            ),
            file=fo
        )
plots(0)

#%% Check which similations failed and why
if False:
    for pname, z in simulations.items():
        bad = 0
        for p in pathways:
            data = pathways[pname].alava[p]
            b, pv = round(data["params"], 2), round(data["p"], 2)
            found = b > 0 and pv < 0.05
            if not found and p in input[pname]:
                print(f'FN: {pname}.{p}: ß={b:5.2f}{"✅" if b > 0 else "❌"}, p={pv:5.2f}{"✅" if pv < 0.05 else "❌"}')
                bad += 1
        if bad:
            print(f"--- {pname}: {bad}/{len(input[pname])} bad")
