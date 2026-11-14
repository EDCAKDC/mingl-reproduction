from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path("results/B004_Ascending")
adata = ad.read_h5ad(
    root / "mingl_neighborhood_results.h5ad"
)

obs = adata.obs.copy()
p = np.asarray(
    adata.obsm["neighborhood_probabilities"]
)

n = len(obs)

assert p.shape[0] == n
assert np.isfinite(p).all()
assert np.allclose(p.sum(axis=1), 1, atol=1e-4)

# MINGL multi-neighborhood candidates
n_positive = (p > 0.25).sum(axis=1)
border = n_positive >= 2

print("Cells:", n)
print("Candidate border cells:", border.sum())
print("Candidate fraction:", border.mean())

# Spatial neighbors: 10 other cells, excluding self
xy = obs[["x", "y"]].to_numpy(dtype=float)
labels = obs["Neighborhood"].astype(str).to_numpy()

if n < 12:
    raise ValueError("At least 12 cells are required")

tree = cKDTree(xy)
_, idx = tree.query(xy, k=12, workers=-1)

# Exclude each query cell even if coordinates are tied
positions = np.broadcast_to(
    np.arange(idx.shape[1]), idx.shape
)
positions = np.where(
    idx == np.arange(n)[:, None],
    idx.shape[1] + 1,
    positions
)

order = np.argsort(positions, axis=1)[:, :10]
neighbors = np.take_along_axis(idx, order, axis=1)

assert np.all(
    neighbors != np.arange(n)[:, None]
)

# Fraction of neighbors with different Neighborhood label
different = (
    labels[neighbors] != labels[:, None]
)

mixing = different.mean(axis=1)

results = pd.DataFrame({
    "x": xy[:, 0],
    "y": xy[:, 1],
    "Neighborhood": labels,
    "n_positive": n_positive,
    "candidate_border": border,
    "mixing_fraction": mixing,
    "max_probability": p.max(axis=1)
})

results.to_csv(
    root / "border_analysis.csv",
    index=False
)

print("\n=== SPATIAL BORDER COMPARISON ===")

for name, mask in [
    ("Single membership", ~border),
    ("Multiple memberships", border)
]:
    values = mixing[mask]

    print(
        f"{name}: n={len(values)}, "
        f"mean_mixing={values.mean():.4f}, "
        f"median_mixing={np.median(values):.4f}"
    )

# Compare border-candidate rates across mixing levels
results["mixing_group"] = pd.cut(
    mixing,
    bins=[-0.001, 0.25, 0.5, 0.75, 1.001],
    labels=["0-25%", "25-50%", "50-75%", "75-100%"]
)

summary = results.groupby(
    "mixing_group", observed=True
).agg(
    n_cells=("candidate_border", "size"),
    candidate_fraction=("candidate_border", "mean")
).reset_index()

print("\n=== BORDER FRACTION BY MIXING ===")
print(summary.to_string(index=False))

summary.to_csv(
    root / "border_mixing_summary.csv",
    index=False
)

# Plot spatial candidate locations
fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(
    xy[~border, 0],
    xy[~border, 1],
    s=1,
    c="lightgray",
    rasterized=True
)

ax.scatter(
    xy[border, 0],
    xy[border, 1],
    s=3,
    c="darkorange",
    rasterized=True
)

ax.set_aspect("equal")
ax.invert_yaxis()
ax.set_title(
    "MINGL multi-neighborhood candidates"
)
ax.set_xlabel("x")
ax.set_ylabel("y")

fig.tight_layout()
fig.savefig(
    root / "candidate_borders.png",
    dpi=300
)
plt.close(fig)

# Plot comparison of spatial mixing
fig, ax = plt.subplots(figsize=(7, 5))

ax.boxplot(
    [mixing[~border], mixing[border]],
    tick_labels=["Single", "Multiple"],
    showfliers=False
)

ax.set_ylabel(
    "Fraction of neighbors with different labels"
)
ax.set_xlabel("MINGL neighborhood memberships")
ax.set_title("Neighborhood mixing comparison")

fig.tight_layout()
fig.savefig(
    root / "border_mixing_comparison.png",
    dpi=300
)
plt.close(fig)

print("\nBORDER ANALYSIS COMPLETE")
