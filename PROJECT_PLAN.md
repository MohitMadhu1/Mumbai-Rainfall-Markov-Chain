# A Markov Chain Model of Monsoon Rainfall in Mumbai

Data-Driven Stochastic Modeling project — final spec.

## Finalised decisions

| # | Decision | Final choice |
|---|----------|--------------|
| 1 | City | Mumbai (IMD grid cell over the city) |
| 2 | States | 3 — Dry / Light / Heavy (2-state Dry/Wet fitted as a baseline) |
| 3 | Season | Main model: monsoon only (Jun–Sep). Comparison: non-monsoon matrix |
| 4 | Thresholds | IMD standard: Dry < 2.5 mm, Light 2.5–64.4 mm, Heavy >= 64.5 mm |
| 5 | Data | IMD gridded daily rainfall 0.25 deg (`pip install imdlib`), 1985–2024. Fit 1985–2014, test 2015–2024 |
| 6 | Extension | Poisson process of Heavy days (>= 64.5 mm) |
| 7 | Validation | All four: order test, simulation check, hold-out forecast, time-homogeneity |
| 9 | Upgrade A | Flood-risk question: P(>= k consecutive Heavy days in a monsoon), return periods, model vs 40 yrs observed |
| 10 | Upgrade B | Monsoon phases: early (Jun), peak (Jul–Aug), late (Sep) matrices; test vs one fixed matrix |
| 11 | Data quality | IMD (government, Pai et al. 2014) cross-checked against NOAA GHCN-Daily gauge for Mumbai |

## Upgrade A — flood risk (headline result)

Mumbai floods when heavy rain persists (26 July 2005). Using the fitted chain:
- Exact probability of a Heavy run of length >= k within a 122-day season
  (Markov-chain recursion / first-passage on an augmented "run length" chain),
  checked against Monte Carlo simulation.
- Expected number of such runs per season and return period T = 1 / P(at least one per season).
- Compare with the observed number of seasons in 1985–2024 containing such a run.

## Upgrade B — monsoon phases

Onset (June) and withdrawal (September) behave differently from the core (Jul–Aug).
- Fit P_early, P_peak, P_late; compare with a single P by likelihood-ratio test and AIC/BIC.
- Compare stationary distributions and mean Heavy-spell lengths across phases.
- Re-run the simulation check with the phase-switching (piecewise time-inhomogeneous) chain.

## Team split (3 people)

| Member | Owns | Report sections |
|--------|------|-----------------|
| A | Data, cleaning, quality cross-check, estimation, monsoon phases (Upgrade B) | Data description, Methodology |
| B | Theory: MLE, P^n, stationary dist., spells, first passage, Poisson, flood-risk maths | Mathematical background, Stochastic model |
| C | Simulator, validation tests 1–4, Poisson tests, flood-risk simulation (Upgrade A) | Simulation results, Validation |
| All | Discussion, limitations, slides, rehearsal | Abstract, Conclusion |
| 8 | Code | Jupyter notebook (runs top to bottom) + `markov.py` core functions, own numpy/scipy code |

## Problem statement

Monsoon rainfall in Mumbai is random day to day, but it is not independent: rainy
days come in spells. Question: can a 3-state Markov chain capture how Mumbai's
monsoon moves between dry, light-rain and heavy-rain days, and how well does it
reproduce spell lengths and heavy-rain frequency seen in 40 years of data?

## Mathematical formulation

- Time: discrete, one step = one day (monsoon days only).
- State space: S = {D, L, H}.
- Process: X_n = state on day n. Assumed time-homogeneous first-order Markov chain.
- Parameters: 3x3 transition matrix P, rows sum to 1 (6 free parameters).
- Poisson part: N(t) = number of Heavy days by day t, rate lambda.

## Theory to derive (report section C)

1. MLE of transition probabilities: p_ij = n_ij / n_i (derive from the likelihood).
2. Confidence intervals for p_ij (asymptotic normal / multinomial).
3. n-step transitions P^n; convergence of P^n rows to pi (eigenvalues of P).
4. Stationary distribution: solve pi P = pi, sum(pi) = 1.
5. Spell length in state i ~ Geometric(1 - p_ii); mean 1/(1 - p_ii).
6. Mean recurrence time m_i = 1/pi_i.
7. Mean first-passage time from D to H (solve linear system).
8. Poisson: counts ~ Poisson(lambda t), gaps ~ Exponential (Geometric in discrete days),
   E[N] = Var[N] = lambda t. Result: counts pass, but gaps fail (Heavy days cluster, p_HH >> pi_H).

## Validation (report section 6)

1. Order test — likelihood-ratio / chi-square: order 0 vs 1 vs 2; AIC/BIC; 2 vs 3 states.
2. Simulation check — 1,000 simulated monsoons: state frequencies vs pi,
   spell-length distributions vs observed (KS / chi-square).
3. Hold-out forecast — one-day-ahead on 2015–2024: Brier score and log-loss vs
   climatology baseline.
4. Time-homogeneity — chi-square test of P(1985–2004) vs P(2005–2024).
5. Poisson — dispersion index test on per-season Heavy counts; gap distribution test.

## Plots

State-transition diagram, transition-matrix heatmap (monsoon vs non-monsoon),
P^n convergence, sample paths (real vs simulated), spell-length histograms
(observed vs geometric vs simulated), hold-out forecast skill, Heavy-day count
plot and inter-arrival histogram, Poisson vs observed count distribution.

## Notebook outline

1. Download + clean IMD data, extract Mumbai cell, classify states
2. Exploratory plots
3. Estimate P (monsoon + non-monsoon), CIs, heatmaps, diagram
4. Theory checks: P^n, pi, spell lengths, recurrence/first-passage times
5. Simulation (own sampler) + sample paths
6. Validation tests 1–4
7. Poisson extension
8. Summary table of results

## Six-week timeline

| Week | Work |
|------|------|
| 1 | Data download, cleaning, state classification, exploratory plots |
| 2 | Estimation, stationary distribution, theory derivations written up |
| 3 | Simulation + validation tests 1–2 |
| 4 | Validation tests 3–4 + Poisson extension |
| 5 | Report draft (10–12 pages) |
| 6 | Slides, rehearse demo, final edits |

## Limitations to discuss

First-order memory only; grid-cell averaging smooths point extremes; no long-term
trends or ENSO/monsoon-onset effects; rainfall amounts not modelled.
Extensions: hidden Markov model, seasonally varying P, gamma model for amounts,
non-homogeneous Poisson trend in Heavy days.

## Key references

- Sonnadara & Jayewardene (2015), Theor. Appl. Climatol., doi:10.1007/s00704-014-1117-z
- Chowdhury et al. (2017), HESS 21, 6541
- Mahanta, Dey & Khosro (2018), Thailand Statistician
- Nandi, Patel & Swain (2024), IMDLIB, doi:10.5281/zenodo.7205414
- Pai et al. (2014), IMD 0.25 deg gridded rainfall dataset, Mausam
