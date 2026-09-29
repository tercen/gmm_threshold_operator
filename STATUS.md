# gmm_threshold_rust_operator — status, 2026-09-29

0.1.0: built from the flowsom_rust_operator skeleton (gather the crosstab, fit, one `.ri`-keyed
table). Gathers the whole crosstab (4 bytes a value), so cohort scale needs the same spill the
other gathering operators need; at 5,000 events per file on the 93-file panel that is 8 M values.
Written for the run-9-on-Tercen workflow (jamie/spectral-pipeline-comparisons, PLAN.md M3).

0.1.1 (2026-09-29): `n_init` restarts (default 10), best log-likelihood kept — written on the assumption that CXCR5's threshold (0.39 vs sklearn's 0.58) was a bad start. It was not: the restarts gave the same 0.39.

0.1.2 (2026-09-29): the real cause was the EM stopping rule. `|Δ mean log-likelihood| < 1e-3` stopped after 6 iterations on CXCR5 (ll −0.5122, means −0.06/0.84); the optimum (ll −0.5100, means −0.05/1.18, threshold 0.56) takes ~30. Tolerance now 1e-6, `max_iter` default 500; every start tried reaches the optimum, so `n_init` is insurance. Fixtures unchanged (the synthetic markers converge either way).
