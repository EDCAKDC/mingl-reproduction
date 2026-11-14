from pathlib import Path
from collections import Counter
import pandas as pd

source = Path.home() / (
    "Downloads/23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv"
)

target_region = "B004_Ascending"
output = Path("data/B004_Ascending_mingl.csv")

cols = [
    "x", "y",
    "unique_region",
    "Cell Type",
    "Neighborhood",
    "Community",
    "Tissue Unit"
]

region_counts = Counter()
parts = []

for chunk in pd.read_csv(
    source,
    usecols=cols,
    chunksize=100000
):
    region_counts.update(
        chunk["unique_region"].value_counts().to_dict()
    )

    selected = chunk[
        chunk["unique_region"] == target_region
    ].copy()

    if not selected.empty:
        parts.append(selected)

print("\n=== LARGEST REGIONS ===")
for name, count in region_counts.most_common(15):
    print(f"{name}: {count:,}")

if not parts:
    raise SystemExit(
        f"Region {target_region} was not found. "
        "Check the region names above."
    )

df = pd.concat(parts, ignore_index=True)

print("\nRegion before filtering:", len(df))
print("\nMissing values:")
print(df.isna().sum())

df = df.dropna(
    subset=["x", "y", "Cell Type", "Neighborhood"]
).reset_index(drop=True)

if len(df) < 300:
    raise SystemExit("Too few valid cells for default MINGL KNN")

df.to_csv(output, index=False)

print("\n=== SELECTED REGION ===")
print("Region:", target_region)
print("Cells:", len(df))
print("Cell types:", df["Cell Type"].nunique())
print("Neighborhoods:", df["Neighborhood"].nunique())

print("\nNeighborhood counts:")
print(df["Neighborhood"].value_counts().to_string())

print("\nSaved:", output.resolve())
print("PREPARATION COMPLETE")
