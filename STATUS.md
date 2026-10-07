# gmm_threshold_rust_operator — status, 2026-09-29

0.1.0: built from the flowsom_rust_operator skeleton (gather the crosstab, fit, one `.ri`-keyed
table). Gathers the whole crosstab (4 bytes a value), so cohort scale needs the same spill the
other gathering operators need; at 5,000 events per file on a 93-file panel that is 8 M values.
Written for a spectral-flow state-analysis workflow on Tercen.

0.1.1 (2026-09-29): `n_init` restarts (default 10), best log-likelihood kept — written on the assumption that one marker's threshold, which differed from sklearn's, came from a bad start. It was not: the restarts gave the same answer.

0.1.2 (2026-09-29): the real cause was the EM stopping rule. `|Δ mean log-likelihood| < 1e-3` stopped after 6 iterations on that marker; the optimum takes ~30. Tolerance now 1e-6, `max_iter` default 500; every start tried reaches the optimum, so `n_init` is insurance. Fixtures regenerated in 0.1.4: the weak bimodal synthetic marker moves by 0.003 with the new start (0.1.3 failed its own install on that).

0.1.3 (2026-09-29): 0.1.2 was the wrong fix. Converging the EM to 1e-6 moved one threshold substantially and flipped another marker to bimodal, neither of which sklearn — and so the reference — does: the early stop at 1e-3 is part of the reference's answer. What the first marker needed was sklearn's **start**: k-means, which in 1-D finds the global 2-means split. The start is now that split, found exactly (`two_means_split`), tolerance back to 1e-3, `max_iter` 100, `n_init` 1. On a 23-marker panel every mode agrees with the sklearn reference and the thresholds agree closely — the same as sklearn itself.

0.1.5 (2026-09-30): memory model sized for this operator — it gathers the crosstab as f32 (4 B/value) plus one f64 subsample copy per marker (≤ 500k × 8 B), so 6 B/value + 60 MB, not FlowSOM's 14 B/value. No code change.
