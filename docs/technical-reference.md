# Scoring and interface reference

## Inputs, normalization and feature selection

`src/run_pipeline.R --help` lists the supported CLI. Run it with `Rscript`.
Metadata are read as text, preserving zero-padded identifiers and additional
columns. Metadata are reordered to matrix-column order; sample sets must match
exactly. Duplicate, missing, blank or control-character identifiers are rejected.
Output paths must be new. Calculations finish before result files are written.

For `--input_scale raw_counts --profile tmm_v1`, counts must be nonnegative
integers with positive libraries. Retain genes with CPM >=1 in at least 20% of
all supplied samples (including follow-up); reset library sizes after filtering; apply edgeR TMM and logCPM with
prior count 1. `tmm_logcpm` may instead be supplied with `tmm_v1`.
VST (variance-stabilizing transformation) sensitivity profile
`vst_sensitivity_v1` requires supplied `deseq2_vst`; the upstream transformation
and design remain the user's responsibility. Normalized inputs must be finite.

Baseline selection requires at least three observations labeled C0 by default.
It retains inclusive type-7 baseline mean quantiles 0.30–0.70, then genes with
baseline MAD >= the eligible genes' type-7 0.60 quantile. MAD uses factor 1.4826.
At least three selected genes are required. The mean/MAD selection step uses baseline only, after all-sample filtering and
normalization. Follow-up counts can therefore influence normalized baseline
values and the resulting selection indirectly.

## Four-component score

For each sample across selected genes:

- IQR uses R's default type-7 quartiles.
- MAD uses the median absolute deviation with factor 1.4826.
- ENT uses `hist(x, breaks=50, plot=FALSE)`, normalizes the returned densities
  to probabilities p, then computes `-sum(p * log(p))` over positive p.
  `breaks` is an adaptive-bin suggestion, not a fixed bin count. Changing
  affine scale can change bin boundaries and therefore entropy.
- Gini first computes `y = scale(x)`, sorts `y - min(y) + 1e-6`, then uses
  `(2*sum(seq_len(n)*y) - (n+1)*sum(y)) / (n*sum(y))`.

A constant sample profile prevents Gini calculation. Nonfinite raw components
or a constant/undefined component across the complete scoring cohort stop
HDI-4 calculation. No imputation is performed.

`z = (raw - mean(raw)) / sd(raw)` uses every scoring sample and sample SD
(denominator n−1). `HDI_4 = (z_IQR + z_MAD + z_ENT + z_Gini)/4`, version
`HDI4-v1`. This CLI does not apply an external training reference. Reference
composition and selected features influence values and ranks.

## Optional BIM and HDI-5

Request `--experimental_bim true`; controls default to `--bim_max_steps 5000
--bim_tolerance 1e-7`. No mixture fitter runs in default mode.
BIM is `1-max(weights)`, a fitted mixture-weight statistic, not evidence of
modality or cell populations. Protocol `returned-poststep-k23-v1` fits k=2 and
k=3 heteroscedastic Gaussian mixtures with equal initial weights, type-7
quantile means spanning 0.15–0.85, initial sample SD and an SD floor of 1e-6.
One step is a full M-step. Convergence requires absolute returned log-likelihood
change below tolerance; the cap is 5,000 steps. Finite returned fits compete on
`(3*k-1)*log(n)-2*logLik`; ties select the first candidate. Only after selection
is the winner rejected for nonconvergence or a final SD at the floor. There is
no alternate-k rescue. Convergence does not establish a global optimum or
stable weights, and synthetic cohorts can legitimately fail.

`HDI_5 = (z_IQR + z_MAD + z_ENT + z_Gini + z_BIM)/5`, version
`HDI5-experimental-v1`, uses the identical four-component reference. Every
sample must have an accepted fit and the complete BIM reference must have
positive finite sample SD. Any rejection sets all HDI-5 to NA with
`unavailable_rejected_fit`; the BIM center/SD are NA. All-accepted constant BIM
retains its mean and zero SD but sets z_BIM and HDI-5 to NA with
`unavailable_constant_BIM`. A raw BIM of zero is valid by itself. No sample
dropping, successful-subset reference, imputation or component omission occurs.
HDI-4 and every row remain intact. Four-only mode uses `not_requested` statuses.

`BIM_fit_records.tsv` records both candidates, returned parameters, criterion,
selection, convergence, iterations, boundary and reasons. Rejected BIM appears
only as `diagnostic_BIM`; the score-table BIM is accepted-only. Reference tables
record centers, SDs, sample counts and ordered IDs. Other products include HVGs,
metadata, roles, retained genes, normalization, configuration, session information,
logs and input/source/output file identities. Generated session information may
contain local environment details; inspect it before sharing outputs.

## R interfaces and synthetic controls

Source modules in this order: `validation.R`, `compute_hvg.R`, `compute_hdi.R`,
`normalize.R`, `bim_prototype.R`, `experimental_scores.R`, `candidate_scores.R`
from `src/`. Prefer `compute_hdi4()` for independent HDI-4 and optional HDI-5.
`freeze_dispersion_reference()` and `score_experimental_dispersion()` support
explicit ordered four/five-component references; HDI-5 requires a named all-TRUE
acceptance vector for the entire reference cohort.

The supported composite API is `compute_hdi4()`. Obsolete `compute_hdi()`,
`compute_hdi_prototype()` and `score_bim_prototype()` scoring wrappers are
removed; shared entropy/Gini helpers, profile validation, returned-parameter
mixture fitting and explicit reference interfaces remain. No supported output
column or version identifier was renamed.

The demo accepts an output parent and always creates a fresh `synthetic_*` run
containing `inputs/` and `results/`. For routine demonstrations use
`Rscript examples/run_demo.R results/demo`; repeat the same command safely.
To explicitly regenerate tracked PNG/PDF figures from the unchanged synthetic
inputs, use `Rscript examples/run_demo.R results/figures --figures`.
Only this explicit option overwrites `figures/synthetic_*.png` and `.pdf`.

`generate_synthetic_data()` supplies the fixed-seed longitudinal demo.
`generate_mixture_scenario()` supplies count-level test controls with known
latent memberships. Counting, normalization and selection can change observed
mixture structure; these controls provide no biological validation.
