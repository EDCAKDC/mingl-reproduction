from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mingl as mg

root = Path("results/B004_Ascending")
outdir = root / "follicle_gradient"
outdir.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(
    root / "mingl_neighborhood_results.h5ad"
)

names = list(map(
    str,
    adata.uns["neighborhood_probability_neighborhoods"]
))
probs = np.asarray(
    adata.obsm["neighborhood_probabilities"],
    dtype=float
)

inner = probs[:, names.index("Inner Follicle")]
outer = probs[:, names.index("Outer Follicle")]

eps = 1e-9

# Same probability-weighted log ratio as MINGL
adata.obs["Score"] = (
    np.log((inner + eps) / (outer + eps))
    * np.maximum(inner, outer)
)

assert np.isfinite(adata.obs["Score"]).all()

print("Computing local gradients for", adata.n_obs, "cells")

result = mg.tl.gb_local_score_gradients(
    adata,
    region_key="unique_region",
    region_value="B004_Ascending",
    score_key="Score",
    k_neighbors=20,
    normalize_by="iqr",
    make_plots=True,
    use_progress=True,
    out_prefix="follicle_local",
    show=False
)

print("\n=== LOCAL GRADIENT SUMMARY ===")
for key, value in result["summary"].items():
    print(key, ":", value)

# Save cell-level gradient measurements
cols = [
    "x", "y", "Neighborhood", "Score",
    "follicle_local_mag",
    "follicle_local_mag_norm",
    "follicle_local_valid"
]
adata.obs[cols].to_csv(
    outdir / "local_gradient_cells.csv",
    index=False
)

# Full tissue map
fig = result["figs"]["scatter"]
ax = fig.axes[0]
ax.invert_yaxis()

fig.savefig(
    outdir / "local_gradient_full.png",
    dpi=300,
    bbox_inches="tight"
)

# Zoom into the prominent lower-right follicle
ax.set_xlim(7000, 8500)
ax.set_ylim(7600, 5700)

fig.savefig(
    outdir / "local_gradient_follicle_zoom.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

# Colorbar
fig = result["figs"]["colorbar"]
fig.savefig(
    outdir / "local_gradient_colorbar.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

print("\nLOCAL GRADIENT ANALYSIS COMPLETE")
