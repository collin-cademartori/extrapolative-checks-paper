# ex2 audit and plan

Findings from reading `ex2/ex2_sim_study.r` and `ex2/ex2_pred_checks.r` against each other and
against ex1, plus the plan for acting on them. Same status key as `TESTING_QUEUE.md`:
**DONE**, **OPEN** (needs a decision), **TEST** (needs a run).

Ordered by how much each item moves a published number, not by when it came up.

---

## 1. Coverage was not measuring coverage — **DONE**

`ex2_sim_study.r` computed its posterior-predictive interval as

```r
y_bounds <- quantile(stat_y_pred, c(0.005, 0.995))          # no indexing
```

`stat_y_pred` is `[draws, T, N]`, so this pooled **every draw, every time point and every unit** into
one global interval — then recomputed that identical interval 240 times inside the `t`/`n` loops.
With unit intercepts spread over `N(4, 3)`, the pooled range is far wider than any single unit's
predictive interval, so `pred_perc` was inflated toward 1 by construction. ex1 has always done this
correctly (`quantile(stat_y_pred[, t, n], ...)`).

The summary's claim that "both models cover in excess of 99%, confirming the check cannot rule out
the intercepts model" was therefore resting on an artifact. The claim may well survive the fix — the
point of Section 5 is that the check *cannot* discriminate — but it has to be re-earned.

**A second bug was hiding behind the first.** The comparison was against `test_ys[t, n]`, the
ORIGINAL column order, while `y_pred` is in the ANCHOR-PERMUTED order the fit used. Only `cor_sq`
and `abs_cors_err` get unpermuted. The global bounds masked the mismatch; fixing the quantile alone
would have introduced a silent unit misalignment. Both fixed together: index per `(t, n)`, compare
against `fit_ys`.

Also added `stopifnot(perm[1] == 1)`. ex2 relies on the treated unit staying in column 1 — for
`delta` and for `unpermute_untreated`'s indexing — without ever checking it. ex1 asserts it.

---

## 2. Absolute error mode — **DONE**

Previously both arms used ratio mode: `err_sd[n] = tau[n] * sigma[n]` with a shared
`tau ~ TN(0.1, 0.05)`. The two arms use **different** `sigma` by design, and that difference leaked
straight into the error scale where nothing justifies it.

Measured over 300 datasets from ex2's own DGP:

```
E[mean_n RMS(y)] = 3.93       E[mean_n sd(y)] = 1.66      ->  2.37x gap
true DGP noise sd = 0.201     (0.1 * max_n sd(latent_n))

  no_ints  sigma = RMS  ->  err prior centre 0.393  =  1.95x the truth
  ints     sigma = sd   ->  err prior centre 0.166  =  0.83x the truth
```

So the **correctly specified arm was handed twice the true observation error** while the
misspecified one got about the right amount — a systematic handicap with no connection to
intercepts, running against the contrast the study exists to measure.

Now: one shared absolute eta for both arms.

```r
sd_y <- apply(fit_ys, 2, sd)
eta_anchor <- mean(sd_y)
ETA_FRAC_EX2 <- 0.12      # 0.201 / 1.66
ETA_CV_EX2   <- 0.5       # same relative width as the old TN(0.1, 0.05)
```

`sd(y_n)` is the anchor because the DGP defines its noise off the latent **sd**, and because it is
arm-neutral — RMS carries the intercepts, which is the leak being closed. `sigma` keeps its
by-design per-arm difference (RMS for `no_ints`, sd for `ints`); eta no longer inherits it.

Absolute mode also collapses eight free per-unit taus to a single shared scalar, which is what
removed ex1's order-statistic funnel and its divergences.

**Added alongside, not separately requested:**

- `eta_med` / `eta_prior_tail` per fit — without them there is no way to see whether the new prior
  is sensible or in conflict, i.e. no way to verify this item landed. Matches ex1.
- `pred_width` — nearly free once per-`(t,n)` bounds exist, and coverage alone cannot distinguish
  "well calibrated" from "intervals so wide they cover everything", which is exactly the failure the
  old global-quantile bug produced. ex1 records it.

---

## 3. Sync `ex2_pred_checks.r` to the study — **TODO**

The prior predictive figure currently describes a model the study does not fit:

| | `ex2_pred_checks.r` | study (`ints` arm) |
|---|---|---|
| `int_scale` | 5 | 3 |
| `int_loc` | absent → 0 | 4 |
| error | `err_scale = 0.2` fixed | estimated |
| `alpha_diag` | absent → 0 | 20 |
| `overall_scales` | `rep(1, N)` | data-driven RMS / sd |

The intercept prior is the substantive one: the figure shows `N(0, 5)`, the study fits `N(4, 3)`.

Note `err_scale = 0.2` in `pred_checks` is almost exactly the DGP truth of 0.201 — the figure had
the right error scale all along and the study did not.

Plan: adopt ex1's named-constant layout (`ETA_FRAC_EX2`, `eta_anchor`, explicit `overall_scales_*`),
work in units where the anchor is 1, and carry the study's multiples unchanged so the two files read
alike and cannot drift.

---

## 4. `ex2_derive_scales.r` — **TODO**

Mirror `ex1/ex1_derive_scales.r` once the constants above are settled: forward-simulate the prior in
plain R, derive `ETA_FRAC_EX2` (and whatever sigma conventions survive), print the derivation, and
assert the committed constants still hold within tolerance so a change to `K_latent`, `T_times`,
`alpha_diag` or a rho prior fails loudly instead of drifting.

---

## 5. Deferred — **OPEN, do not act without discussion**

- **A `no_ints` prior predictive panel.** `ex2_pred_checks.r` only shows the `ints` arm, though the
  paper's contrast is between the two.
- **The DGP's noise anchor.** `0.1 * max_n sd(latent_n)` uses `max`, the least stable anchor
  available — the same instability that steered ex1 toward mean-based anchors. Measured, the true
  noise sd ranges **0.103 to 0.446** across datasets, so the "truth" itself varies four-fold rep to
  rep. Everything downstream inherits that as extra variance.
- **`f_treat_sd <- 1.9` is hardcoded** as the post-treatment divergence of `f_alt`, while the
  realised `arima.sim(ar = 0.9, n = 30)` sd varies around it. The treatment-window separation is
  therefore a fixed absolute amount rather than one tied to each dataset's realised scale.
- **99% vs 95% intervals.** ex2 reports 99%, ex1 reports 95%. ex2's label matches its code; the
  question is only whether the two examples should agree.

---

## Remaining ex1/ex2 divergences not covered above

Recorded so they are not lost, none of them urgent:

- `anchor_order(test_ys, K_latent)` in ex2 vs `K_latent + 1` in ex1. ex1 takes the extra column as
  deliberate slack; ex2 does not.
- ex2 has no per-unit overfitting measure analogous to ex1's `noise_abs_tr`. May simply not apply —
  ex2's story is about location/correlation entanglement, not about absorbing noise.

---

## Related: `TESTING_QUEUE.md` A1 is now resolved

A1 asked whether ex1's DGP should set `alpha_diag` to match its fits. It now does — the ex1 DGP
passes `alpha_diag = 20` and `absolute_error = TRUE`, so the DGP and the nonstationary arm share
their loadings prior, their rho prior (both `Beta(8, 2)`) and their error parametrization. That
entry can be closed.
