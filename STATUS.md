# gmm_threshold_rust_operator — status, 2026-09-29

0.1.0: built from the flowsom_rust_operator skeleton (gather the crosstab, fit, one `.ri`-keyed
table). Gathers the whole crosstab (4 bytes a value), so cohort scale needs the same spill the
other gathering operators need; at 5,000 events per file on the 93-file panel that is 8 M values.
Written for the run-9-on-Tercen workflow (jamie/spectral-pipeline-comparisons, PLAN.md M3).

0.1.1 (2026-09-29): `n_init` restarts (default 10), best log-likelihood kept — written on the assumption that CXCR5's threshold (0.39 vs sklearn's 0.58) was a bad start. It was not: the restarts gave the same 0.39.

0.1.2 (2026-09-29): the real cause was the EM stopping rule. `|Δ mean log-likelihood| < 1e-3` stopped after 6 iterations on CXCR5 (ll −0.5122, means −0.06/0.84); the optimum (ll −0.5100, means −0.05/1.18, threshold 0.56) takes ~30. Tolerance now 1e-6, `max_iter` default 500; every start tried reaches the optimum, so `n_init` is insurance. Fixtures regenerated in 0.1.4: the weak bimodal synthetic marker moves by 0.003 with the new start (0.1.3 failed its own install on that).

0.1.3 (2026-09-29): 0.1.2 was the wrong fix. Converging the EM to 1e-6 moved IgM's threshold from 1.06 to 1.81 and flipped NKP44 to bimodal, neither of which sklearn — and so the reference — does: the early stop at 1e-3 is part of the reference's answer. What CXCR5 needed was sklearn's **start**: k-means, which in 1-D finds the global 2-means split. The start is now that split, found exactly (`two_means_split`), tolerance back to 1e-3, `max_iter` 100, `n_init` 1. On the Lyme panel: 23/23 modes agree with run 9, median |diff| 0.004, none above 0.1 — the same as sklearn itself.
