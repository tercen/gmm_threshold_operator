# gmm_threshold_rust_operator

One positivity threshold per marker, from a two-component Gaussian mixture. The rule is the
CYTOSHRINK state analysis's (`state_analysis.py`, section 14): fit the mixture to a seeded
subsample of cells; if the two means are at least `min_separation` apart the marker is
**bimodal** and the threshold is the midpoint of the means; otherwise it is **unimodal** and the
threshold is `mean_low + tail_sd × sd_low`.

| | |
|---|---|
| rows | markers — one threshold per row |
| columns | cells |
| y | the transformed value (asinh, normalised) |
| output | one row per marker, joined on the row: `threshold`, `mode`, `separation`, `mean_low`, `sd_low`, `mean_high`, `sd_high`, `weight_low`, `pct_positive` (over all cells), `cells_fitted`, `em_iterations` |
| image | `ghcr.io/tercen/gmm_threshold_rust_operator` |

## Properties

| property | default | meaning |
|---|---|---|
| `subsample` | 500000 | Cells drawn, seeded and without replacement, to fit each marker; the same cells for every marker. 0 uses all. |
| `seed` | 5024 | Seed of that draw. |
| `min_separation` | 0.5 | Means at least this far apart: bimodal. |
| `tail_sd` | 2.0 | Unimodal threshold: `mean_low + tail_sd × sd_low`. |
| `max_iter` | 100 | EM iterations at most. |
| `n_init` | 1 | EM starts per marker, best log-likelihood kept (sklearn's `n_init`). Extra starts at seeded quantile pairs. Not needed for parity: the default start is the exact 1-D 2-means split, which is what sklearn's k-means seed reaches. |

## Parity

The mixture is fitted by EM from the exact 1-D 2-means split (what sklearn's k-means seed reaches), to
a per-sample log-likelihood tolerance of 1e-3 with a variance floor of 1e-6 — sklearn's
`GaussianMixture(2, n_init=1)` rule, early stop included: converging further changes the answer on
weakly bimodal markers and the reference did not. On a 23-marker spectral panel the
thresholds match the sklearn reference closely with every mode agreeing. `cargo test` checks a clear bimodal mixture, a unimodal one, and the component ordering.

Feed the thresholds back as a row factor to a fraction-above step to get percent positive per
group; see `fraction_above_rust_operator`.
