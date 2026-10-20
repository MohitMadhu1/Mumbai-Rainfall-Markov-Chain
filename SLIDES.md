# Presentation plan — 15 minutes (11 talk + 4 Q&A)

One idea per slide. Every number below is in `analysis.ipynb` / `report_draft.md`.
Figures are in `figures/`. Suggested speaker split in brackets (A / B / C).

| # | Time | Slide | Show | Say (one line) |
|---|------|-------|------|----------------|
| 1 | 0:30 | **Title** — A Markov Chain Model of Monsoon Rainfall in Mumbai | team names | — |
| 2 | 0:30 | **The problem** [A] | `02_example_season.png` | Rain comes in spells, not independent days. How long do they last, and how often does heavy rain persist? |
| 3 | 1:00 | **Data** [A] | `01_data_quality.png` | 40 years of official IMD gridded rainfall, no missing days. A NOAA rain gauge agrees: seasonal r = 0.94. |
| 4 | 1:00 | **The model** [B] | `04_transition_diagram.png` | 3 states from IMD thresholds. Tomorrow depends only on today. 6 numbers describe the monsoon. |
| 5 | 1:00 | **Theory** [B] | `05_n_step_convergence.png` + small table | π = (0.37, 0.56, 0.07) matches the data. The chain forgets its start within a week. Heavy spells last 1.5 days on average. |
| 6 | 1:00 | **Methodology and code** [C] | **live notebook**: the estimation cell, then the simulation cell | MLE p̂ = n_ij / n_i, our own sampler, 1,000 synthetic monsoons. |
| 7 | 1:00 | **Does it fit? (1)** [C] | `07_spell_lengths.png` | Heavy spells fit (p = 0.79); dry spells do not, because of monsoon breaks. |
| 8 | 1:00 | **Does it fit? (2)** [C] | `09_day_of_season_phases.png` + forecast table | 15.5% better than climatology on 2015–2024. Onset, peak and withdrawal matrices fix the seasonal shape. |
| 9 | 1:00 | **Headline: flood risk** [A] | `10_flood_risk.png` | Exact answer: 54% chance of a 3-day Heavy run per monsoon (observed 55%). It rose from 0.44 to 0.63 after 2005. |
| 10 | 0:45 | **Poisson twist** [B] | `11_poisson.png` (middle panel) | Counts look Poisson, but 39% of gaps are one day vs 8% expected. Heavy rain clusters, and the Markov chain captures it. |
| 11 | 0:45 | **Conclusions and limitations** [all] | 3 bullets + 3 bullets | Works for persistence and risk. Misses dry-spell memory, amounts and the climate trend. |
| 12 | 1:00 | **Individual contributions** | table | A: data, estimation, phases · B: theory, Poisson, flood maths · C: simulation, validation, forecasts |
| — | 4:00 | **Q&A** | keep `analysis.html` open | — |

## Likely questions and short answers

- **Why 3 states and not 2?** Rows Light and Heavy of P differ (LRT p ≈ 10⁻³⁶). After Heavy, P(Heavy) = 0.35; after Light, 0.08.
- **Is first order enough?** It gives 90% of the likelihood gain. Second order is significant but improves forecasts by only 2.6% with 3× the parameters.
- **Why IMD and not a single station?** It is official, peer-reviewed (Pai et al. 2014) and has no gaps. The gauge has 9% missing days and was used to cross-check.
- **Why did dry spells fail?** Monsoon "breaks": dry spells are a mix of 1-day pauses and long breaks, and a geometric spell cannot be both.
- **Is the increase real?** Yes, on three counts: the time-homogeneity test (p = 0.003); seasons with a 3-day run went from 7/20 to 15/20 (Fisher's exact p = 0.025); and the independent Santacruz gauge shows the same rise, more strongly (Heavy days 7.4% → 11.2%, trend p = 0.001), so it is not an artefact of the IMD grid. Consistent with Goswami et al. (2006) and Roxy et al. (2017).
- **Isn't your forecast skill just seasonality?** No. Against a baseline that knows the day of the season, the chain still has about 10% Brier skill.
- **Why no KS test?** Spell lengths and counts are whole numbers; KS p-values are invalid for discrete data, so we use binned chi-square.
- **Isn't the 54% vs 55% match in-sample?** Yes. So we also fitted on 1985–2014 only and predicted 2015–2024: 3-day runs expected in 5.2 of 10 seasons, observed 8; 4-day runs 2.3 vs 2; 5-day 0.9 vs 1. No binomial test rejects (all p > 0.1). The 3-day under-prediction matches the rising-risk finding.
- **Is it a good forecaster?** No, and it isn't meant to be: 15.5% skill over climatology comes from persistence alone. Its value is describing patterns and risk.
- **What about 26 July 2005?** It was one extreme day (≈ 490 mm). Our states do not separate 65 mm from 490 mm, so it is a limitation and an extension (amount model).
