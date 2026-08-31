# Open questions and pending tests

Everything raised and deferred while getting the sampler under control. Ordered within each
section by how much it matters, not by when it came up. Merge sequencing lives separately in
`MERGE_PLAN.md`; this file is about what still needs deciding or measuring.

Status key: **OPEN** (needs a decision), **TEST** (needs a run), **WATCH** (data will arrive on its
own from the next production run).

---

## A. Model specification — changes the posterior, belongs in the paper

### A1. `alpha_diag` in ex1's prior-predictive DGP — **OPEN**

Every *fit* in both studies now uses `alpha_diag = 20` (zero-avoiding inverse-gamma on the loading
diagonal). ex1's DGP never sets it, so it generates under `sample_model`'s default of `0`, the plain
half-normal.

So the nonstationary arm — the one the example treats as correctly specified — is not fitting under
the prior its own data was generated from, for the loading diagonal. This predates today (nonstat
was at 10, DGP at 0); syncing the fits to 20 only made it conspicuous.

Two coherent positions:

- **Leave the DGP at 0.** The inverse-gamma is a computational device to repel collapsed-loading
  modes, not a claim about the world; "correctly specified" refers to structure. Costs nothing,
  regenerates no data. *Current state, and my weak preference.*
- **Set the DGP to 20.** nonstat then genuinely fits under the generating prior. But this
  regenerates every dataset and invalidates comparison with every run so far.

Either way the choice needs one sentence in the paper. Does not arise in ex2, whose DGP is a
separate adversarial process that never touches `ife_named`'s priors.

### A2. Re-separate ex1's `stat_weak` and `stat_strong` — **OPEN, blocks the merge**

The two stationary arms now produce nearly the same results, so the example no longer demonstrates
the contrast it exists to show. Their only difference is the tau prior: `N(0.1, 0.1)` vs
`N(0.05, 0.05)`. Expected to be a small change, but it is a *production* change to
`ex1_sim_study.r`, so the merge should wait for it (see `MERGE_PLAN.md`, precondition 2).

### A3. `K_latent` sensitivity — **TEST**

Current state is *not* uniform:

| | DGP K | fitted K |
|---|---|---|
| ex1 nonstat | 4 | 4 |
| ex1 stat_weak / stat_strong | 4 | **5** (`K_latent + 1`) |
| ex2 both | 3 | 3 |

Left as-is deliberately. The argument for the asymmetry: ex1's stationary fits have no intercept and
no factor means, so the level must come from the factors and an extra factor gives it somewhere to
go — the same structural argument that justifies ex1's 2x RMS sigma multiple. Both ex2 fits have an
explicit level mechanism, so an extra factor there adds an unidentified direction and nothing else.

Wanted as a sensitivity analysis regardless. The clean version fits K-1 / K / K+1 with the DGP held
fixed, as `ex1_k_sweep.r` already does, and reports `rhat_M`, `cor_sq` and the delta error.

### A4. Anchor ordering into ex1 — **OPEN**

ex2 permutes columns before fitting (`anchor_order`) so the leading K x K loading block is full rank;
ex1 does not. Worth adding for consistency, but not a one-liner: ex1 indexes `noise_abs_tr`,
`y_means` and the treated unit by column position, so it needs `unpermute_untreated`-style mapping
back. It touches the metric feeding the overfitting plot, so it wants care.

---

## B. Sampler configuration — computation only, decide empirically

### B1. `max_treedepth` in ex2 — **WATCH**

Dropped from 12 to ex1's default of 10 in `aa54c4c`. The 4-task smoke then showed the `ints` arm
hitting the limit on up to 5% of transitions (`n_tree` up to 483). Not yet evidence: that smoke ran
at `iter = 150 / warm = 150`, and poor adaptation inflates saturation on its own.

If saturation persists at production warmup, the 12 was compensating for something real in the
intercepts geometry and should go back. `n_tree` is recorded per fit, so this is a lookup.

### B2. Are the escalation thresholds sane for ex2? — **WATCH**

`rhat_M > 1.01`, `ess_delta1 < 400` and a 0.1% divergence rate were all calibrated on ex1 and
inherited unchanged. ex2 differs in ways that could matter: intercepts/factor means, 4 chains
instead of 3, 500 sampling iterations instead of 2000.

Raw draw count is *not* the right lens here — ESS is driven by autocorrelation, which may differ
between the studies in either direction. The question is empirical: what fraction of ex2's 1600 fits
escalate, and on which criterion. Worth checking after the first ~100 tasks before committing to the
full grid, since escalation cost scales with it.

### B3. `n_chains`: 3 (ex1) vs 4 (ex2) — **OPEN**

Arbitrary. ex1 dropped to 3 to speed up computation. Plan is to put ex1 back to 4 later. Affects
R-hat quality directly, so worth doing before the runs that get reported.

### B4. Base warmup length — **TEST**

ex1 warms up 500 and samples 2000; ex2 does 500/500. There is no theoretical ratio requirement —
warmup and sampling have different jobs and each needs only to be long enough for its own. But
warmup was *measured* to improve `rhat_M` monotonically in all three `stat_strong` cases at fixed
sampling length (1.092 -> 1.038, 1.017 -> 1.004, 1.014 -> 1.009), with `ess_delta1` moving in step.

So the open question is not the ratio but whether ex1's *base* warmup of 500 is simply too short.
Raising it may be cheaper than refitting escalated fits at 8000. One caveat already seen: 265
stat_strong picked up 8 divergences at `warm = 4000` having had none at 500 or 2000.

### B5. Divergence escalation threshold — **TEST**

`ESCALATE_DIV_RATE = 0.001` means >6 divergences at base. In the last run 37% of `stat_weak` fits
had at least one divergence but only 11% crossed the threshold, so roughly a quarter of fits keep
1-6 divergences that are never escalated away. Those draws are biased in principle.

Lowering the threshold spends compute only where divergences appear (unlike raising base
`adapt_delta`, which taxes every fit), but would move `stat_weak` from ~11% to ~37% escalating.
Worth checking first whether those low-count fits differ in delta error — `n_div` is now recorded
per fit, so this is answerable from the next run without any new sampling.

### B6. Fresh vs fixed seed across escalation rounds — **OPEN**

Each escalated round draws a new seed, so a refit is an independent re-roll rather than a
continuation. That is what made the old intermediate 4000-iteration rung useless (it improved
`rhat_M` in 12 of 21 trajectories and worsened it in 6). Removing the marginal rung mitigated it, but
the underlying choice was never revisited. A fixed seed would make escalation monotone in
computation; it would also correlate the rounds.

---

## C. Diagnostics and presentation

### C1. Does low E-BFMI matter? — **TEST, answerable from recorded data**

E-BFMI is recorded but deliberately not an escalation criterion: it is inert to sampling length by
construction, was measured to be nearly inert to warmup, and across ex1's stationary round-1 fits is
essentially unrelated to estimand mixing (`cor(ebfmi, rhat_M) = -0.12`, `cor(ebfmi, ess_delta1) =
+0.07`). It is also strikingly seed-dependent: unit 278 stat_weak scored 0.19 in one run and 0.52
under an identical configuration differing only in seed.

The check that settles whether it can be ignored: compare delta error and `noise_abs_tr` between the
low-E-BFMI fits (~7% of stationary fits, `ebfmi_min < 0.3`) and the rest. `ebfmi_min` is recorded
per fit, so no new sampling is needed.

### C2. The `rho` <-> `||Lambda||` hypothesis for low E-BFMI — **TEST**

The untested structural explanation. In ex1 sigma is fixed and `omega_sq = 1`, so the remaining
multiplicative freedom is between a near-unit-root factor's realised amplitude and the loading
magnitude, under a steep `beta(97,3)` prior. Consistent with nonstat being clean (flatter
`beta(8,2)`, fits on differences). A fixed-rho run would test it in about ten minutes. Only worth
doing if C1 says E-BFMI actually tracks something.

### C3. Reporting the escalation ladder honestly — **OPEN**

A fit stops escalating the moment it passes, so the *recorded* `rhat_M` distribution is
stopped-on-success and reads low. The 70-rep rebuild showed this concretely: every completed rep
ended at `rhat_M <= 1.010`, 0 of 70 above threshold in all three arms — which is a real result but
not an unbiased sample of the diagnostic.

The `n_rounds` distribution must be reported alongside it. Already noted in the shared code comment
in `sample_model.r`; needs to make it into the paper.

### C4. Sensitivity to elevated-R-hat fits — **OPEN**

The original plan before escalation existed: check how the summaries (delta error, overfitting
metric) move as fits with elevated `rhat_M` are included or held out. Now a secondary robustness
check rather than the main defence, since escalation resolves most of them — but still the right
answer for whatever survives a full ladder.

### C5. SBC for delta — **OPEN, probably out of scope**

Raised as a way to make the "R-hat on loadings is benign" argument rigorous rather than
argumentative. Simulation-based calibration on delta would show the estimand is calibrated
regardless of loading-configuration mixing. Expensive, and `rhat_M` plus the conditional
independence argument may be enough. Parked.

---

## D. Verification runs pending

### D1. `ex1_memory_check.r` on the desktop — **TEST**

Verifies the OOM mechanism and, more usefully, sets the worker count. Laptop numbers (heap metric):
`iter=1000` bare 254 MB vs named 197 MB; `iter=2500` bare 457 vs named 271. Post-fix cost scaled as
~148 MB + 50 KB per sampling iteration, i.e. ~736 MB per fit at `iter = 12000`.

```
cd ex1 && Rscript ex1_memory_check.r 2000,8000,12000 <n_workers>
```

On Linux the RSS column reads `VmHWM` and is trustworthy; on macOS ignore it and read the heap
column. Falsification criteria are in the script header.

### D2. First real ex2 run — **TEST**

Never run under the ported escalation ladder, `alpha_diag = 20`, or `max_treedepth = 10`. Only a
4-task reduced smoke has exercised the plumbing. Feeds B1 and B2.

### D3. `origin/stan-tweak` — **OPEN**

Remote-only branch, no local checkout, never examined. Determine whether anything in it still
matters before the merge.

---

## E. Bookmarked 2026-08-30: does extra K widen or narrow the overfitting gap?

**TEST**, after we understand behaviour at the true K.

Fitting more latent factors than the truth lowers the rank-k residual floor, which is the quantity
coverage tracks (coverage holds when resid_sd / err_sd ~ 1). Measured over 200 datasets, true K = 4,
true noise sd 2.0:

| fitted k | best rank-k residual sd | vs rank-4 |
|---|---|---|
| 3 (K−1) | 1.312 | 126% |
| 4 (true) | 1.045 | — |
| 5 (K+1) | 0.801 | 77% |
| 6 (K+2) | 0.567 | 54% |
| 7 (K+3) | 0.340 | 33% |

So K is a **second dial**, orthogonal to `err_sd`: extra factors create slack that lets `err_sd` go
lower — more overfitting, higher S1 p — while keeping predictive coverage near nominal. k=6 would
support roughly the `err_sd = 0.5` configuration with coverage intact instead of 0.775.

**The open question:** if K rises for the STATIONARY arms, should it also rise for `nonstat`, and
what happens to the overfitting gap between them?

Collin's hypothesis: `nonstat` has a mitigation the stationary model lacks. It models the
*increments*, so its fitted trajectory is an integrated process and far smoother than an AR(1) path
at the level. Both arms overfit more as K grows, but the stationary arm should overfit more, and the
gap may widen because the noise being absorbed is rough while the true factors are smooth.

A sharper version of the same mechanism, which makes the "widens" prediction more likely: absorbing
iid LEVEL noise requires strongly NEGATIVELY autocorrelated increments (alternating up-down). The
nonstationary arm puts an AR(0.8) prior on its increments, which pushes toward POSITIVE
autocorrelation — so its extra capacity is spent against a prior that actively resists the shape the
absorption requires. The stationary arm faces no such obstruction: its factors need only level-scale
wiggle, and AR(0.97) with unit marginal variance permits steps of ~0.24 in the factor, which times
sigma ~8 is ~2 at the data scale — comparable to the noise sd of 2. So extra factors should buy the
stationary arm much more absorption per factor than they buy `nonstat`.

**Test:** sweep k in {4, 5, 6, 7} for both arms, holding `err_sd` fixed, and record `noise_abs_tr`,
`pred_perc`, S1 p and the sampler diagnostics. Two things to watch: whether the `noise_abs_tr` gap
widens or narrows, and whether `rhat_loadings` degrades — the K vs K+1 comparison over 22 paired
datasets found no detectable difference, but K+2 and K+3 are extrapolation from that.

If the gap widens, K becomes the more natural dial to lead with in the paper: fitting more factors
than the truth is an ordinary modelling choice, whereas a hand-set error scale needs justifying.
