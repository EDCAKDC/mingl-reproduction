#!/usr/bin/env bash
set -euo pipefail

# Run from the extracted publication package directory.
# Uses only results that already exist in the user's local MINGL directory.
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UPSTREAM="${MINGL_SOURCE_DIR:-$HOME/projects/Mingl}"
REPOSITORY="${GITHUB_REPOSITORY:-EDCAKDC/mingl-reproduction}"

printf '\n=== Assemble independent MINGL reproduction repository ===\n'
printf 'Source: %s\nOutput: %s\nRepository: %s\n' "$UPSTREAM" "$BASE" "$REPOSITORY"

if [[ ! -d "$UPSTREAM" || ! -d "$UPSTREAM/src/mingl" ]]; then
  echo 'ERROR: Could not find the official MINGL source directory.' >&2
  echo 'Set MINGL_SOURCE_DIR=/path/to/Mingl and rerun.' >&2
  exit 1
fi

if ! command -v gh >/dev/null 2>&1; then
  echo 'ERROR: GitHub CLI (gh) is missing.' >&2
  exit 1
fi
if ! gh auth status >/dev/null 2>&1; then
  echo 'ERROR: GitHub CLI authentication is not configured.' >&2
  echo 'Run: gh auth login --hostname github.com --git-protocol https --web' >&2
  exit 1
fi

mkdir -p "$BASE/scripts" "$BASE/figures" "$BASE/tables" "$BASE/env"

# Fail early instead of publishing a repository that refers to missing scripts.
scripts=(
  prepare_mingl_region.py
  run_mingl_b004.py
  analyze_mingl_borders.py
  run_mingl_border_enrichment.py
  run_mingl_follicle_gradient.py
  run_mingl_local_gradient.py
  fix_mingl_local_gradient.py
)
for name in "${scripts[@]}"; do
  if [[ ! -f "$UPSTREAM/$name" ]]; then
    echo "ERROR: Required script not found: $UPSTREAM/$name" >&2
    exit 1
  fi
done
for name in "${scripts[@]}"; do
  cp "$UPSTREAM/$name" "$BASE/scripts/$name"
  echo "Copied script: $name"
done

# Copy only selected, publication-suitable PNGs. Avoid the invalid old
# IQR-normalized local gradient images, raw data and per-cell outputs.
copy_asset() {
  local relative_source="$1"
  local destination="$2"
  local source="$UPSTREAM/results/B004_Ascending/$relative_source"
  if [[ -f "$source" ]]; then
    cp "$source" "$BASE/$destination"
    echo "Copied: $destination"
  else
    echo "Not found (skipped): $source" >&2
  fi
}

copy_asset 'neighborhood_memberships.png' 'figures/neighborhood_memberships.png'
copy_asset 'candidate_borders.png' 'figures/candidate_borders.png'
copy_asset 'border_mixing_comparison.png' 'figures/border_mixing_comparison.png'
copy_asset 'border_celltype_enrichment.png' 'figures/border_celltype_enrichment.png'
copy_asset 'border_count_legend.png' 'figures/border_count_legend.png'
copy_asset 'border_celltype_legend.png' 'figures/border_celltype_legend.png'
copy_asset 'follicle_gradient/follicle_spatial_clusters.png' 'figures/follicle_spatial_clusters.png'
copy_asset 'follicle_gradient/follicle_cell_composition.png' 'figures/follicle_cell_composition.png'
copy_asset 'follicle_gradient/follicle_gradient_scores.png' 'figures/follicle_gradient_scores.png'
copy_asset 'follicle_gradient/local_gradient_corrected_full.png' 'figures/local_gradient_corrected_full.png'
copy_asset 'follicle_gradient/local_gradient_corrected_zoom.png' 'figures/local_gradient_corrected_zoom.png'

copy_asset 'border_mixing_summary.csv' 'tables/border_mixing_summary.csv'
copy_asset 'border_pair_counts.csv' 'tables/border_pair_counts.csv'
copy_asset 'follicle_gradient/cluster_summary.csv' 'tables/official_gradient_cluster_summary.csv'
copy_asset 'follicle_gradient/cluster_probability_bin_fc.csv' 'tables/cluster_probability_bin_fc.csv'
copy_asset 'follicle_gradient/steepness.json' 'tables/steepness.json'

# Record exact upstream commit without publishing the official code itself.
if git -C "$UPSTREAM" rev-parse --verify HEAD >/dev/null 2>&1; then
  git -C "$UPSTREAM" rev-parse HEAD > "$BASE/env/upstream_mingl_commit.txt"
fi

# Required figures referred to directly in README must be present.
required_figures=(
  candidate_borders.png
  border_mixing_comparison.png
  follicle_spatial_clusters.png
  follicle_cell_composition.png
  local_gradient_corrected_zoom.png
)
for name in "${required_figures[@]}"; do
  if [[ ! -s "$BASE/figures/$name" ]]; then
    echo "ERROR: README figure missing: $BASE/figures/$name" >&2
    echo 'The package was assembled but NOT published. Restore the missing PNG and rerun.' >&2
    exit 1
  fi
done

cd "$BASE"
if [[ ! -d .git ]]; then
  git init -b main
fi

git add .gitignore README.md docs scripts figures tables env publish.sh
if ! git diff --cached --quiet; then
  git commit -m 'Add independent MINGL reproduction, figures and analysis notes'
fi

# The user explicitly runs this script; creating a public repo is the intended action.
if gh repo view "$REPOSITORY" >/dev/null 2>&1; then
  echo 'GitHub repository already exists; pushing to main.'
  if ! git remote get-url origin >/dev/null 2>&1; then
    git remote add origin "https://github.com/${REPOSITORY}.git"
  fi
  git push -u origin main
else
  echo 'Creating a new PUBLIC GitHub repository and pushing main.'
  gh repo create "$REPOSITORY" --public --source=. --remote=origin --push
fi

printf '\nPUBLISHED: https://github.com/%s\n' "$REPOSITORY"
git status --short --branch
