# Figure notes

- **candidate_borders.png** — Cells with two or more MINGL membership probabilities above 0.25. These are candidates for interfaces, not independently verified biological borders.
- **border_mixing_comparison.png** — Comparison of original neighborhood-label disagreement among the ten nearest *other* cells. This is an internal consistency check since the original labels also contributed to the membership model.
- **border_celltype_enrichment.png** — Official MINGL enrichment figure for the Outer/Inner Follicle pair, generated with `plot_border_enrichment` and a 0.25 membership threshold. Counts and relative compositions depend on the original annotations.
- **follicle_spatial_clusters.png** — Five transition clusters (3,522 cells). Colors correspond to ranks based on mean probability score, not raw cluster IDs. The plot was generated with `use_quantiles=False`.
- **follicle_cell_composition.png** — B and CD4+ T cell fractions across ranked clusters. A composition transition is evident, though CD4 fraction is not strictly monotonic.
- **local_gradient_corrected_full.png** — Local gradient magnitude, range-normalized for all cells but colored only where either Inner or Outer Follicle probability exceeds 0.25.
- **local_gradient_corrected_zoom.png** — Enlarged view of a prominent follicle. The raw local regression was computed on the entire region; selection was only for display. The color upper bound is the 99th percentile among displayed cells.

The earlier `local_gradient_full.png` and `local_gradient_follicle_zoom.png` used unstable IQR normalization (effective denominator 10⁻⁸ for an almost-zero full-section score IQR) and were deliberately excluded.
