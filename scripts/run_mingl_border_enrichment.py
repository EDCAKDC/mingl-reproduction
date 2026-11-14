from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mingl as mg

out = Path("results/B004_Ascending")
adata = ad.read_h5ad(
    out / "mingl_neighborhood_results.h5ad"
)

p = np.asarray(
    adata.obsm["neighborhood_probabilities"]
)

names = [
    str(x) for x in
    adata.uns["neighborhood_probability_neighborhoods"]
]

assert p.shape[1] == len(names)
assert "Cell Type" in adata.obs

positive = p > 0.25
cell_types = adata.obs["Cell Type"].astype(str).to_numpy()

rows = []

for i, j in combinations(range(len(names)), 2):
    a = positive[:, i]
    b = positive[:, j]

    both = a & b
    left = a & ~b
    right = b & ~a

    if both.sum() == 0:
        continue

    # Cell types represented in all three subsets
    shared = 0
    for ct in np.unique(cell_types):
        mask = cell_types == ct
        if (
            np.sum(mask & both) >= 5
            and np.sum(mask & left) >= 5
            and np.sum(mask & right) >= 5
        ):
            shared += 1

    rows.append({
        "neighborhood_1": names[i],
        "neighborhood_2": names[j],
        "both": int(both.sum()),
        "only_1": int(left.sum()),
        "only_2": int(right.sum()),
        "shared_cell_types": shared
    })

pairs = pd.DataFrame(rows).sort_values(
    "both", ascending=False
)

pairs.to_csv(
    out / "border_pair_counts.csv",
    index=False
)

print("\n=== TOP NEIGHBORHOOD PAIRS ===")
print(pairs.head(15).to_string(index=False))

eligible = pairs[
    (pairs.both >= 30)
    & (pairs.only_1 >= 50)
    & (pairs.only_2 >= 50)
    & (pairs.shared_cell_types >= 1)
]

if eligible.empty:
    raise SystemExit(
        "No pair has enough cell-type support for "
        "the default enrichment plot. "
        "Inspect border_pair_counts.csv."
    )

top = eligible.iloc[0]
n1 = top["neighborhood_1"]
n2 = top["neighborhood_2"]

print("\n=== SELECTED BORDER ===")
print("Neighborhood 1:", n1)
print("Neighborhood 2:", n2)
print("Joint membership cells:", top["both"])

figures = mg.pl.plot_border_enrichment(
    adata=adata,
    n1=n1,
    n2=n2,
    cell_type_col="Cell Type",
    pos_threshold=0.25,
    min_count=5,
    show=False,
    label_dots=True,
    label_fontsize=7
)

for fig, name in zip(
    figures,
    [
        "border_celltype_enrichment.png",
        "border_count_legend.png",
        "border_celltype_legend.png"
    ]
):
    fig.savefig(
        out / name,
        dpi=250,
        bbox_inches="tight"
    )
    plt.close(fig)

print("\nMINGL BORDER ENRICHMENT COMPLETE")
