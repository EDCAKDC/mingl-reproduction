from pathlib import Path
import numpy as np
import pandas as pd
import mingl as mg

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path("results/B004_Ascending")
out.mkdir(parents=True, exist_ok=True)

# 1. Load one complete spatial region
cells = mg.pp.read_file(
    "data/B004_Ascending_mingl.csv"
)

print("Cells:", cells.n_obs)
print("Neighborhoods:",
      cells.obs["Neighborhood"].nunique())

assert cells.n_obs >= 300
assert cells.obs["Neighborhood"].nunique() >= 2

# 2. Calculate neighborhood composition centroids
print("\n=== CALCULATING CENTROIDS ===", flush=True)

centroids = mg.tl.centroid_Calculation(
    cells,
    k=10,
    cluster_col="Cell Type",
    neighborhood_col="Neighborhood",
    region_col="unique_region"
)

print("Centroid shape:", centroids.shape, flush=True)

# 3. Compute neighborhood membership probabilities
print("\n=== CALCULATING PROBABILITIES ===", flush=True)

mg.tl.cpu_gmm_probability(
    CELLS_ADATA=cells,
    CENTROIDS_ADATA=centroids,
    cluster_col="Cell Type",
    neighborhood_col="Neighborhood",
    region_key="unique_region",
    k=10,
    num_processes=4
)

probs = np.asarray(
    cells.obsm["neighborhood_probabilities"]
)

print("Probability matrix:", probs.shape)

row_sums = probs.sum(axis=1)
print(
    "Rows with zero total probability:",
    int((row_sums == 0).sum())
)

print(
    "Rows with nonfinite probabilities:",
    int((~np.isfinite(probs)).any(axis=1).sum())
)

# Stop rather than interpret invalid probabilities
if not np.isfinite(probs).all() or (row_sums == 0).any():
    raise RuntimeError(
        "Invalid or zero probability rows. "
        "Inspect GMM results before continuing."
    )

# 4. Identify cells with multiple positive memberships
print("\n=== IDENTIFYING BORDER CANDIDATES ===", flush=True)

mg.tl.findPositives(
    cells,
    prob_key="neighborhood_probabilities",
    threshold=0.25,
    result_key="Count_Above_Threshold"
)

counts = cells.obs[
    "Count_Above_Threshold"
].value_counts().sort_index()

print("\nPositive neighborhood counts:")
print(counts.to_string())

# 5. Save the analytical results
cells.write_h5ad(
    out / "mingl_neighborhood_results.h5ad",
    compression="gzip"
)

# 6. Plot the spatial distribution
obs = cells.obs

fig, ax = plt.subplots(figsize=(9, 8))

scatter = ax.scatter(
    obs["x"],
    obs["y"],
    c=obs["Count_Above_Threshold"],
    s=1,
    cmap="viridis",
    rasterized=True
)

ax.set_aspect("equal")
ax.invert_yaxis()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(
    "MINGL: Positive neighborhood memberships"
)

fig.colorbar(
    scatter,
    ax=ax,
    label="Neighborhood count (P > 0.25)"
)

fig.tight_layout()
fig.savefig(
    out / "neighborhood_memberships.png",
    dpi=250
)

plt.close(fig)

print("\nSaved outputs to:", out.resolve())
print("MINGL INITIAL ANALYSIS COMPLETE")
