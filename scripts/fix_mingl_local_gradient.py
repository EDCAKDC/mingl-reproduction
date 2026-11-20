from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path("results/B004_Ascending")
out = root / "follicle_gradient"

df = pd.read_csv(out / "local_gradient_cells.csv")
adata = ad.read_h5ad(
    root / "mingl_neighborhood_results.h5ad"
)

xy = df[["x", "y"]].to_numpy(dtype=float)
original_xy = adata.obs[["x", "y"]].to_numpy(dtype=float)

assert len(df) == adata.n_obs
assert np.allclose(xy, original_xy)

names = list(map(
    str,
    adata.uns["neighborhood_probability_neighborhoods"]
))
p = np.asarray(
    adata.obsm["neighborhood_probabilities"],
    dtype=float
)

inner = p[:, names.index("Inner Follicle")]
outer = p[:, names.index("Outer Follicle")]

score = df["Score"].to_numpy(dtype=float)
raw_gradient = df[
    "follicle_local_mag"
].to_numpy(dtype=float)

q25, q75 = np.percentile(score, [25, 75])
score_iqr = q75 - q25
score_range = np.ptp(score)

# Same scaling used by MINGL's normalize_by="range"
nn_dist, _ = cKDTree(xy).query(
    xy, k=2, workers=-1
)
median_nn = np.median(nn_dist[:, 1])

if not np.isfinite(score).all() or not np.isfinite(raw_gradient).all():
    raise ValueError("Non-finite score or raw-gradient values in input CSV")
if not np.isfinite(p).all():
    raise ValueError("Non-finite neighborhood probabilities in AnnData")
if not np.isfinite(score_range) or score_range <= 0:
    raise ValueError("Score range must be finite and positive for normalization")
if not np.isfinite(median_nn) or median_nn <= 0:
    raise ValueError("Median nearest-neighbor distance must be finite and positive")
norm_gradient = (
    raw_gradient * median_nn / score_range
)

# Only highlight cells relevant to the two follicles
# 0.25 matches our earlier membership threshold.
support = np.maximum(inner, outer) > 0.25
valid = support & np.isfinite(norm_gradient)

assert valid.sum() > 0

print("\n=== SCORE DIAGNOSTICS ===")
print("Score Q25:", q25)
print("Score Q75:", q75)
print("Actual Score IQR:", score_iqr)
print("Score range:", score_range)
print("Exact-zero Score fraction:", np.mean(score == 0))
print(
    "Inner/Outer total probability < 1e-8:",
    np.mean((inner + outer) < 1e-8)
)

print("\n=== CORRECTED GRADIENT ===")
print("Follicle-relevant cells:", int(valid.sum()))
print("Median:", np.median(norm_gradient[valid]))
print("Mean:", np.mean(norm_gradient[valid]))
print("P95:", np.percentile(norm_gradient[valid], 95))
print("P99:", np.percentile(norm_gradient[valid], 99))

# Use a robust display limit; numerical values are not changed
vmax = np.percentile(norm_gradient[valid], 99)
vmax = max(vmax, 1e-12)

fig, ax = plt.subplots(figsize=(11, 9))

ax.scatter(
    xy[:, 0], xy[:, 1],
    s=0.6, c="lightgray",
    alpha=0.55, rasterized=True
)

sc = ax.scatter(
    xy[valid, 0],
    xy[valid, 1],
    c=norm_gradient[valid],
    s=5,
    cmap="inferno",
    vmin=0,
    vmax=vmax,
    rasterized=True
)

ax.set_aspect("equal")
ax.invert_yaxis()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(
    "Inner/Outer Follicle local gradient\n"
    "Range-normalized; follicle-relevant cells"
)

cb = fig.colorbar(sc, ax=ax)
cb.set_label("Range-normalized local gradient")

fig.tight_layout()
fig.savefig(
    out / "local_gradient_corrected_full.png",
    dpi=300
)

ax.set_xlim(6900, 8500)
ax.set_ylim(7600, 5500)

fig.savefig(
    out / "local_gradient_corrected_zoom.png",
    dpi=300
)

plt.close(fig)

df["range_normalized_gradient"] = norm_gradient
df["follicle_relevant"] = support
df.to_csv(
    out / "local_gradient_corrected.csv",
    index=False
)

print("\nCORRECTED LOCAL GRADIENT COMPLETE")
