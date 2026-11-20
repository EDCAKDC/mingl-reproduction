# MINGL analysis of a human colon region

**Zhixuan Zheng · Independent reproduction notes · November 2026**

## Research question

Spatial transcriptomics and multiplexed imaging provide each cell's identity and position, but conventional neighborhood labels place cells into discrete categories. MINGL describes tissue organization in terms of continuous neighborhood membership probabilities. I used one public CODEX colon region to understand how those probabilities relate to mixed neighborhoods and to gradients across lymphoid follicles.

This work focused on reproducing and interpreting selected MINGL analysis steps. It was not designed as a new-method benchmark, a re-annotation of the original tissue, or a test of clinical significance.

## Data and software

The input was `23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv`, a public merged CODEX intestine table with coordinates and annotations. I selected `unique_region == "B004_Ascending"`, retaining 21,232 cells. The metadata used here were `x`, `y`, `unique_region`, `Cell Type`, `Neighborhood`, `Community`, and `Tissue Unit`.

Analysis used the Hickey Lab MINGL Python package, version 0.1.0, installed from its public repository on Python 3.11 / Pop!_OS. Each cell's original `Neighborhood` label was treated as a supplied annotation. Neither cell types nor neighborhood identities were inferred independently during this project.

## Step 1: Neighborhood membership probabilities

I used the official `centroid_Calculation` and `cpu_gmm_probability` functions. The first summarizes the cell-type composition of nearby cells for the existing neighborhood categories. The second compares each cell's local composition against those reference profiles to assign a probability vector across neighborhoods. A neighborhood size of `k=10` was used.

The resulting probability matrix had shape **21,232 × 19**. No row contained non-finite probability values or a zero probability sum.

Using the official `findPositives` function and threshold `P > 0.25`, 19,295 cells had one above-threshold membership; 1,900 had two; 37 had three. Therefore, **1,937 / 21,232 = 9.12%** had two or more above-threshold memberships. This is a probabilistic boundary-candidate definition, not a histological ground truth.

## Step 2: Spatial mixing consistency check

I independently computed each cell's nearest ten *other* cells from the supplied XY coordinates. For each cell, I measured the fraction of neighbors assigned a different original `Neighborhood` label. The mean was:

| Membership category | Cells | Mean different-neighborhood neighbor fraction | Median |
|---|---:|---:|---:|
| One positive neighborhood | 19,295 | 0.3161 | 0.30 |
| Two or more positive neighborhoods | 1,937 | 0.5220 | 0.50 |

The fraction of MINGL multi-membership candidates increased with the amount of local mixing:

| Fraction of nearest neighbors with different labels | Cells | Multi-membership candidate rate |
|---|---:|---:|
| 0–25% | 9,702 | 3.659% |
| 25–50% | 6,072 | 10.672% |
| 50–75% | 3,047 | 15.852% |
| 75–100% | 2,411 | 18.706% |

The pattern is compatible with MINGL's intended interpretation, but the two analyses share the original neighborhood labels. The check therefore tests internal spatial coherence rather than external biological validity.

**Threshold sensitivity (code-only follow-up):** The candidate definition depends on the `P > 0.25` cutoff. The analysis script now prints candidate counts at 0.20, 0.25, and 0.30 as a descriptive check. Those extra counts have **not been run or reported as new results here**. An external boundary reference and tests across donors would still be required for validation.

## Step 3: Pairwise neighborhood interfaces

I enumerated pairs of neighborhood probabilities that simultaneously exceeded 0.25. The five highest pair counts were:

| Pair | Cells with both memberships |
|---|---:|
| Outer Follicle / Inner Follicle | 174 |
| Transit Amplifying Zone / Secretory Epithelial | 173 |
| Transit Amplifying Zone / Plasma Cell Enriched | 155 |
| Plasma Cell Enriched / CD8+ T Enriched IEL | 98 |
| Stroma & Innate Immune / Stroma | 92 |

The official `plot_border_enrichment` method ran successfully for the Outer/Inner Follicle pair. A cell can contribute to more than one neighborhood-pair count. Pair counts were not corrected for neighborhood prevalence or spatial interface length.

## Step 4: Inner/Outer Follicle gradient

The region included 1,240 cells labeled Inner Follicle and 2,282 labeled Outer Follicle. I called MINGL's `mingl_neighborhoods_scverse` with the two follicle labels, `k=20`, five probability bins, `n_clusters=5`, `distance_max=1000`, `use_quantiles=False`, and `random_state=42`.

The probability-weighted score was calculated as

\[
S_i = \log \left(\frac{P_i(\mathrm{Inner})+\varepsilon}{P_i(\mathrm{Outer})+\varepsilon}\right) \max\{P_i(\mathrm{Inner}),P_i(\mathrm{Outer})\}.
\]

The five clusters were displayed from higher to lower mean score. These displayed ranks are *not* the raw K-means cluster identifiers.

| Rank, Inner-like → Outer-like | n | Mean score | B-cell fraction | CD4+ T-cell fraction |
|---|---:|---:|---:|---:|
| 0 | 611 | +19.7013 | 0.8936 | 0.0753 |
| 1 | 660 | +5.0043 | 0.4682 | 0.3212 |
| 2 | 752 | −4.8759 | 0.1569 | 0.4668 |
| 3 | 784 | −7.8531 | 0.1314 | 0.6645 |
| 4 | 715 | −14.5050 | 0.1399 | 0.5203 |

In the spatial plots, the ordered clusters formed visually interpretable layers around prominent follicles. The Inner-like cells had a high B-cell fraction, whereas the outer clusters showed a higher contribution from CD4+ T cells. The CD4+ fraction was not strictly monotonic. Although the gradient ranking follows the score by construction, the changes in annotated cell composition were observed independently of that ordering rule.

Using MINGL's official `steepness_score` on the saved cell-level results gave **steepness = 5.1315589948** and **normalized slope = 0.4220005396**. The official ordered original cluster labels were `['1','4','3','0','2']`, with mean scores `[19.7013, 5.0043, −4.8759, −7.8531, −14.5050]`.

This steepness metric is based on differences between ordered cluster mean scores, not uniform physical-distance intervals. Its absolute value depends on analysis settings and should not be compared to a different figure or sample without matching the analysis protocol.

## Step 5: Local gradient diagnostic

The official `gb_local_score_gradients` calculated a spatially local linear regression of the Inner/Outer Score using `k_neighbors=20` for all 21,232 cells. All cells had valid local fits.

Its default IQR normalization was **not interpretable for this sample**:

| Diagnostic | Value |
|---|---:|
| Actual full-tissue Score IQR | 6.0772 × 10⁻¹⁷ |
| Fraction of exact-zero Scores | 49.20% |
| Fraction with P(Inner)+P(Outer) < 10⁻⁸ | 76.40% |
| Score range | 41.4465 |
| Median nearest-neighbor distance | 15.2643 |

MINGL imposes a lower IQR bound of 10⁻⁸, so its default normalized values became extremely large for some cells (mean approximately 49.4 million). This was a scaling problem, not evidence of an extraordinarily strong biological gradient.

I retained the original locally fitted raw gradient magnitudes but recalculated the normalization using `raw_gradient × median_nearest_neighbor_distance / full_score_range`. For visualization, I then displayed only cells with `max(P(Inner), P(Outer)) > 0.25`. The fit still used all local neighbors; the restriction was applied **only to displayed cells**.

This selection retained **3,455 cells**, with range-normalized gradient median **0.04713**, mean **0.05365**, 95th percentile **0.13134**, and 99th percentile **0.17066**. Higher values appeared around parts of the follicular interface. This image is an exploratory diagnostic and **not an exact reproduction of the original paper panel or default normalization**.

## Interpretation

The observed spatial layering, cell-composition transition, and neighborhood membership overlap are consistent with the biological intuition behind MINGL. The project also highlighted a practical failure mode of global IQR normalization when the target pair of neighborhoods is rare relative to the full tissue.

It remains a **one-region, one-donor** case study. It does not show that MINGL outperforms other tools, demonstrate new cell states, or independently validate interfaces. The neighborhood identities were previously supplied, and the exact publication figures may use different tissues and parameters. Future work should reproduce the published protocol on matching source data, compare independent regions/donors, and examine parameter sensitivity before drawing stronger conclusions.

## Files and provenance

- `scripts/`: analysis commands used locally, separate from the official source.
- `tables/`: compact aggregated results and diagnostic summary; no per-cell information.
- `figures/`: selected results, excluding the invalid IQR-normalized plots.
- `docs/reproducibility.md`: how to obtain the original package/data and run the analysis.

**Reference:** Van Batavia K et al. *MINGL: Quantifies Borders, Gradients, and Heterogeneity in Multicellular Tissue Organization*. bioRxiv (2026). https://www.biorxiv.org/content/10.64898/2026.03.24.713296v1
