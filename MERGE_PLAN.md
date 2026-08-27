# Merge plan: parametrization-diagnostics → sampling-fixes → main

Status: agreed in principle, **not started**. Nothing below has been executed. Blocked until
(a) the ex1 sampler run completes cleanly and (b) ex1's weak/strong prior separation is settled —
see "Preconditions".

## Branch roles

- **`main`** — public-facing reproduction code only.
- **`sampling-fixes`** — the staging bridge. Diagnostic and experimental branches land here first;
  only the curated result goes on to `main`. Keeps the full working record.
- **`parametrization-diagnostics`** — current work. +14 over `sampling-fixes`.
- **`dgp-diag`** — the original sampler investigation (where Pathfinder init came from). Forked
  from `main`, never rebased. Its keepable content was *manually* ported into `sampling-fixes`,
  which is why git records no ancestry link. Leave as-is.
- **`unitnorm-estsigma`** — unit-norm loadings + estimated sigma. Dead end, superseded by the
  revert to `ife_named` in `32738ed`. Leave as-is; ignore for merging.
- **`origin/stan-tweak`** — remote-only, no local branch. Unexamined.

## Decision: curate onto sampling-fixes (option 3)

`parametrization-diagnostics` is **74 files, +8471 / −125** over `sampling-fixes`, but splits
sharply by purpose:

```
production   :  8 files,   +583 / -125    <- the real work
scaffolding  : 66 files,  +7888           <- 9 Stan variants, 57 probe scripts
```

93% of the branch is exploration testing hypotheses that were mostly refuted. That record has real
value — several probe scripts are the evidence behind comments now in the production code — but it
belongs on `sampling-fixes`, not on `main`.

Plan: build a clean commit series on `sampling-fixes` covering only the eight production files, and
keep `parametrization-diagnostics` as an archived reference branch (do not delete).

### The eight production files

| file | what it carries |
|---|---|
| `sample_model.r` | CSV cleanup; split rhat (`rhat_M`/estimands/loadings/`cor_sq`); `escalation_ladder()` + `fit_with_escalation()` |
| `pathfinder_init.r` | Pathfinder CSV cleanup |
| `ex1/ex1_sim_study.r` | 2× RMS sigma anchor; `noise_abs_tr`; `EX1_LADDER`; checkpointing; `stat_strong` truncated-normal tau |
| `ex1/ex1_sim_study_summary.r` | overfit plot keyed on `noise_abs_tr` |
| `ex1/ex1_pred_checks.r` | synced with the study's config |
| `ex1/ex1_rsq_priors.r` | synced with the study's config |
| `ex2/ex2_sim_study.r` | `EX2_LADDER`; per-fit diagnostics; checkpointing |
| `.gitignore` | checkpoint dirs, `progress.log` |

## What NOT to do

**Do not rebase onto `99b0375` and drop the later `sampling-fixes` commits.** This was considered
and ruled out on evidence. Survival of each commit's added lines at
`parametrization-diagnostics` HEAD:

| commit | file | survival |
|---|---|---|
| `99b0375` | `pathfinder_init.r` | 42/42 (100%) |
| | `sample_model.r` | 46/49 (93%) |
| | `ife_named.stan` | 15/15 (100%) |
| `0ba81bb` | `ex1_sim_study.r` | 4/10 (40%) |
| `0756610` | `ex2_sim_study.r` | 60/63 (95%) |
| `92b2bfd` | summary bugfix | 3/3 (100%) |
| `ce9bb59` | summary bugfix | 2/2 (100%) |

`0756610` carries the ex2 DGP rewrite (`T_times` 20→30, AR 0.96→0.9, `N_unc`, `K_unc`, noise
0.05→0.1, `sim` parameterised, the `seq_len` zero-group fix, the `Y` orientation fix) and the
Cholesky anchor ordering (`anchor_order`, `unpermute_untreated`).
**`parametrization-diagnostics` does not touch any of it** — zero matching diff lines — so dropping
`0756610` would lose work nothing reintroduces.

`0ba81bb`'s 40% looks alarming but is not: five of the six non-surviving lines are cosmetic
(`exp_vars` was extended; `n_chains`/`seed` moved into the ladder's argument list). Exactly one
semantic reversal: `err_scale = 0.1`, the Dirac-delta `stat_strong` prior, replaced by the
truncated normal in `3df7512`. Its load-bearing content (`source("../pathfinder_init.r")`,
`alpha_diag`, `fit_scales = FALSE`, `pathfinder_init = TRUE`) is all still live.

## Correction to the working narrative

The recollection that "the original fix had ex1 and ex2 share a value of sigma" is not quite what
the history shows, and the real footprint is smaller:

- Fixing sigma was **never reverted**. `fit_scales = FALSE`/`0` is still in place at HEAD (3 sites
  in ex1, 2 in ex2).
- The two examples never shared a literal value: ex1 used `sd(y)`; ex2 used `RMS(y)` for `no_ints`
  and `sd(y)` for `ints`. What they shared was the *convention* — anchor sigma to a plain,
  unmultiplied data-scale statistic.
- What actually changed is **one line in ex1**:
  `apply(test_ys, 2, sd)` → `2 * apply(test_ys, 2, function(y) sqrt(mean(y^2)))`, justified by the
  realised-amplitude correction for near-unit-root factors in fits that have no level mechanism.
  ex2 keeps 1× because both its fits have one.

## Preconditions before executing

1. ex1 long run completes without crashing, and the `n_rounds` / `rhat_M` distributions look
   acceptable.
2. ex1's `stat_weak` vs `stat_strong` priors are re-separated — they currently produce nearly the
   same results, so another production commit to `ex1_sim_study.r` is expected. Merge after that
   lands, not before.
3. ex2 verified under the ported escalation ladder and the inherited thresholds (`rhat_M > 1.01`,
   `ess < 400`, div rate > 0.1%, ad floor 0.95), which were calibrated on ex1's sampler and may not
   suit ex2's intercepts / `alpha_diag = 0` / 4 chains.

## Open questions

- Does any scaffolding deserve promotion to `sampling-fixes` as a permanent record — e.g. the
  prior-predictive check behind the 2× RMS multiple — rather than living only in branch history?
- What is `origin/stan-tweak`, and does anything in it still matter?
