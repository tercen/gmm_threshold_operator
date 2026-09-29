# gmm_threshold_rust_operator — status, 2026-09-29

0.1.0: built from the flowsom_rust_operator skeleton (gather the crosstab, fit, one `.ri`-keyed
table). Gathers the whole crosstab (4 bytes a value), so cohort scale needs the same spill the
other gathering operators need; at 5,000 events per file on the 93-file panel that is 8 M values.
Written for the run-9-on-Tercen workflow (jamie/spectral-pipeline-comparisons, PLAN.md M3).
