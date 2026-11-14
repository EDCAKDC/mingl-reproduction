from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mingl as mg

root = Path("results/B004_Ascending")
out = root / "follicle_gradient"
out.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(
    root / "mingl_neighborhood_results.h5ad"
)

names = [
    str(x) for x in
    adata.uns["neighborhood_probability_neighborhoods"]
]

p = np.asarray(
    adata.obsm["neighborhood_probabilities"]
)

inner = "Inner Follicle"
outer = "Outer Follicle"

assert inner in names and outer in names
assert np.isfinite(p).all()

# Official gradient function expects probability columns
# to be present in adata.obs.
adata.obs[inner] = p[:, names.index(inner)]
adata.obs[outer] = p[:, names.index(outer)]

print("Cells:", adata.n_obs)
print("Inner Follicle cells:",
      int((adata.obs["Neighborhood"] == inner).sum()))
print("Outer Follicle cells:",
      int((adata.obs["Neighborhood"] == outer).sum()))

print("\n=== RUNNING MINGL GRADIENT ===", flush=True)

adata, df_sub, fc, windows, edges = (
    mg.tl.mingl_neighborhoods_scverse(
        adata,
        tu1=inner,
        tu2=outer,
        target_neighborhoods=(inner, outer),
        ks=(20,),
        k=20,
        distance_max=1000,
        n_bins=5,
        use_quantiles=False,
        n_clusters=5,
        random_state=42
    )
)

# Preserve numeric cluster ids, but order their
# presentation by average Inner-to-Outer score.
obs = adata.obs.copy()
sub = obs[
    obs["Neighborhood"].isin([inner, outer])
].copy()

sub["cluster"] = sub[
    "Probability_Bin_Cluster"
].astype(str)

assert not sub["cluster"].isna().any()
assert len(sub) > 0

order = (
    sub.groupby("cluster")["Score"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

rank = {c: i for i, c in enumerate(order)}
sub["gradient_rank"] = sub["cluster"].map(rank)

sub["is_B"] = (
    sub["Cell Type"].astype(str) == "B"
).astype(int)

sub["is_CD4T"] = (
    sub["Cell Type"].astype(str) == "CD4+ T cell"
).astype(int)

summary = (
    sub.groupby("gradient_rank")
    .agg(
        n_cells=("Score", "size"),
        mean_score=("Score", "mean"),
        B_fraction=("is_B", "mean"),
        CD4T_fraction=("is_CD4T", "mean")
    )
    .reset_index()
    .sort_values("gradient_rank")
)

summary.to_csv(
    out / "cluster_summary.csv",
    index=False
)
sub[
    [
        "x", "y", "Neighborhood", "Cell Type",
        "Score", "Probability_Level",
        "cluster", "gradient_rank"
    ]
].to_csv(
    out / "gradient_cells.csv",
    index=False
)
fc.to_csv(out / "cluster_probability_bin_fc.csv")

print("\n=== FOLLICLE GRADIENT SUMMARY ===")
print(summary.round(4).to_string(index=False))

# Figure 1: spatial arrangement of transition clusters
fig, ax = plt.subplots(figsize=(9, 8))

ax.scatter(
    obs["x"], obs["y"],
    s=0.4, color="lightgray",
    rasterized=True
)

sc = ax.scatter(
    sub["x"], sub["y"],
    c=sub["gradient_rank"],
    cmap="viridis",
    vmin=0, vmax=4,
    s=5, rasterized=True
)

ax.invert_yaxis()
ax.set_aspect("equal")
ax.set_title("Inner-to-Outer Follicle transition clusters")
ax.set_xlabel("x")
ax.set_ylabel("y")

cb = fig.colorbar(sc, ax=ax, ticks=range(5))
cb.set_label("Cluster rank (Inner-like to Outer-like)")

fig.tight_layout()
fig.savefig(out / "follicle_spatial_clusters.png", dpi=300)
plt.close(fig)

# Figure 2: B-cell and CD4 T-cell compositions
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(
    summary["gradient_rank"],
    summary["B_fraction"],
    marker="o", label="B cells"
)
ax.plot(
    summary["gradient_rank"],
    summary["CD4T_fraction"],
    marker="o", label="CD4+ T cells"
)

ax.set_xlabel("Transition cluster rank")
ax.set_ylabel("Fraction of cells")
ax.set_xticks(range(5))
ax.set_title("Cell composition across follicle gradient")
ax.legend()

fig.tight_layout()
fig.savefig(out / "follicle_cell_composition.png", dpi=300)
plt.close(fig)

# Figure 3: probability-weighted log ratio
fig, ax = plt.subplots(figsize=(8, 5))

data = [
    sub.loc[sub.gradient_rank == i, "Score"].to_numpy()
    for i in range(5)
]

ax.boxplot(data, showfliers=False)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(range(5))
ax.axhline(0, linestyle="--", color="gray")
ax.set_xlabel("Transition cluster rank")
ax.set_ylabel("Probability-weighted log ratio")
ax.set_title("Inner / Outer Follicle score distribution")

fig.tight_layout()
fig.savefig(out / "follicle_gradient_scores.png", dpi=300)
plt.close(fig)

print("\nGRADIENT ANALYSIS COMPLETE")
print("Results:", out.resolve())
