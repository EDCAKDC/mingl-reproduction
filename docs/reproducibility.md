# Reproducibility

## Inputs

The original public merged CODEX intestine dataset is expected at:

`~/Downloads/23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv`

It contains per-cell protein measurements, coordinates and annotations. Only these columns are used in the MINGL analysis:

`x`, `y`, `unique_region`, `Cell Type`, `Neighborhood`, `Community`, `Tissue Unit`.

The analyses select only the region `B004_Ascending` (21,232 cells in the tested dataset). The raw merged table and the extracted region-level CSV are not distributed here. See the MINGL paper and [original source repository](https://github.com/HickeyLab/Mingl) for dataset attribution and download instructions.

## Environment

The analysis used Python 3.11 and MINGL 0.1.0. For a new installation:

```bash
conda create -n mingl python=3.11 -y
conda activate mingl

git clone https://github.com/HickeyLab/Mingl.git ~/projects/Mingl
python -m pip install -e ~/projects/Mingl
python -c 'import mingl; print(mingl.__version__)'
```

For exact reproduction, compare the Git revision recorded in `env/upstream_mingl_commit.txt` (if present). Because upstream code may change, the installed version alone may not capture all implementation differences.

## Run

Run the following from the root of **this independent repository**, after downloading the dataset and installing MINGL:

```bash
mkdir -p data results/B004_Ascending
python scripts/prepare_mingl_region.py
python scripts/run_mingl_b004.py
python scripts/analyze_mingl_borders.py
python scripts/run_mingl_border_enrichment.py
python scripts/run_mingl_follicle_gradient.py
python scripts/run_mingl_local_gradient.py
python scripts/fix_mingl_local_gradient.py
```

The first script reads the large merged CSV in chunks and writes a much smaller local region CSV. The following scripts run the analysis. The local `results/` directory is ignored by Git; selected figures and aggregate summaries are copied into `figures/` and `tables/` for publication.

The spatial mixing check uses **10 other cells**, excluding the cell itself. The MINGL centroid/GMM routines use the installed package's `k=10` semantics. The gradient analysis uses `k=20`, five clusters, `random_state=42`, and `use_quantiles=False`. Those choices matter when comparing results against different versions of the manuscript.

### Checks and interpretation

- `analyze_mingl_borders.py` prints **descriptive** candidate-border counts at probability cutoffs 0.20, 0.25, and 0.30. All saved figures and originally reported counts still use the 0.25 cutoff. This sensitivity output is a newly added check, **not a rerun or independent validation**.
- `fix_mingl_local_gradient.py` checks that scores, raw gradients, probabilities, score range, and median neighbor distance are finite and that normalization denominators are positive. Its alternative display scales the previously fitted raw gradient by `median_neighbor_distance / (max(Score) - min(Score))`.
- The plotting filter `max(P_inner, P_outer) > 0.25` affects **displayed points only**, not which cells were included in the local regression. The output is an exploratory diagnostic, not a replacement for the official MINGL score or published figure.
- To verify new analytical results after these documentation and script improvements, **rerun the scripts with the original data locally**. No new dataset-level measurements are claimed from the code-only updates.

## What is and is not in this GitHub project

Included: analysis scripts written for this reproduction, a summary of actual results, small aggregate CSVs, and selected PNG figures.

Excluded: the 2.91 GB merged source CSV, extracted region CSV, the 21,232-cell AnnData output, per-cell probability matrices, cell-level gradient CSV files, console logs, and raw model intermediates. This avoids duplicating the published data and keeps the repository easy to inspect.

No new GPU training was performed during these MINGL neighborhood and follicle analyses; they used CPU implementations. This project is separate from the STELLAR GraphSAGE GPU reproduction.
