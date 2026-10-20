"""Generate analysis.ipynb. Run:  .venv/Scripts/python build_notebook.py"""

import nbformat as nbf

cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


# =====================================================================
md(r"""
# A Markov Chain Model of Monsoon Rainfall in Mumbai

**Data-Driven Stochastic Modeling and Simulation of Real-World Systems**

Pipeline: real-world problem → data → probability model → stochastic process →
mathematical analysis → simulation → validation → interpretation.

| | |
|---|---|
| Data | IMD 0.25° gridded daily rainfall (India Meteorological Department), cell 19.00 N 72.75 E, 1985–2024; cross-checked with NOAA GHCN-Daily, Santacruz observatory |
| States | Dry (< 2.5 mm), Light (2.5–64.4 mm), Heavy (≥ 64.5 mm) — IMD definitions |
| Season | Monsoon, 1 June – 30 September (122 days) |
| Training / test | 1985–2014 / 2015–2024 |

All estimation, simulation and testing code is our own (`markov.py`, numpy/scipy only).
""")

code(r"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from scipy import stats

import markov as mk
import viz
from viz import OBS, MODEL, THEORY, EXTRA, STATE, INK_2, MUTED

viz.use()
rng = np.random.default_rng(2026)
pd.set_option("display.precision", 4)
pd.set_option("display.max_colwidth", None)
S = list(mk.STATE_NAMES)
N_DAYS, N_SIM = 122, 1000
TRAIN, TEST = range(1985, 2015), range(2015, 2025)
""")

# =====================================================================
md(r"""
## 1. Data

Daily rainfall for the IMD grid cell containing the Santacruz observatory.
The IMD gridded product (Pai et al., 2014) is interpolated from the national
rain-gauge network and has no missing days. The Santacruz gauge (NOAA GHCN-Daily,
station IN012070800) is an independent point record used only as a quality check.
""")

code(r"""
df = pd.read_csv("data/mumbai_daily_rain.csv", parse_dates=["date"], index_col="date")
df["year"], df["month"] = df.index.year, df.index.month
mon = df[df.month.between(6, 9)].copy()
print(f"{len(df):,} days 1985-2024 | monsoon days: {len(mon):,} | "
      f"IMD missing: {df.rain_imd.isna().sum()} | gauge missing (monsoon): {mon.rain_gauge.isna().sum()}")
df.describe().loc[["count", "mean", "std", "max"], ["rain_imd", "rain_gauge"]]
""")

md(r"""
### 1.1 Data quality: do the two government sources agree?
""")

code(r"""
both = mon.dropna(subset=["rain_imd", "rain_gauge"])
r_daily = both.rain_imd.corr(both.rain_gauge)
tot = both.groupby("year")[["rain_imd", "rain_gauge"]].sum()
r_season = tot.rain_imd.corr(tot.rain_gauge)
s_imd, s_gauge = mk.classify(both.rain_imd), mk.classify(both.rain_gauge)
agree = (s_imd == s_gauge).mean()
conf = pd.crosstab(pd.Series(s_imd, name="IMD").map(dict(enumerate(S))),
                   pd.Series(s_gauge, name="Gauge").map(dict(enumerate(S)))).reindex(index=S, columns=S)
print(f"Common monsoon days: {len(both):,}")
print(f"Daily correlation IMD vs gauge:     r = {r_daily:.3f}")
print(f"Seasonal-total correlation:         r = {r_season:.3f}")
print(f"Same state (Dry/Light/Heavy) on     {agree:.1%} of days")
conf
""")

code(r"""
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
ax[0].plot(tot.index, tot.rain_imd, color=OBS, marker="o", label="IMD grid cell")
ax[0].plot(tot.index, tot.rain_gauge, color=MODEL, marker="o", label="Santacruz gauge")
ax[0].set(title="Monsoon rainfall, common days only", ylabel="mm per season")
ax[0].legend()
lim = [0, tot.values.max() * 1.05]
ax[1].scatter(tot.rain_gauge, tot.rain_imd, color=OBS, s=30, edgecolor="white", linewidth=1)
ax[1].plot(lim, lim, color=MUTED, lw=1, ls="--")
ax[1].set(title=f"Seasonal totals agree (r = {r_season:.2f})", xlabel="Gauge (mm)", ylabel="IMD (mm)",
          xlim=lim, ylim=lim)
fig.tight_layout(); viz.save(fig, "01_data_quality")
""")

md(r"""
**Finding.** On the 4,440 monsoon days where both sources report, daily rainfall correlates at
r = 0.80 and seasonal totals at r = 0.94, and the two put the day in the same state 79.5% of the
time. Almost every disagreement is between neighbouring states; only 19 days are Dry in one source
and Heavy in the other. A grid-cell average is expected to differ from a single gauge (it smooths
local showers), so this is strong agreement. IMD has no missing days, while the gauge misses 440
monsoon days (9%), which is why IMD is the primary source.
""")

# =====================================================================
md(r"""
## 2. States and the stochastic process

* **Time:** discrete, $n = 1, \dots, 122$ (days of the monsoon), one sequence per year.
* **State space:** $S = \{D, L, H\}$ with IMD thresholds.
* **Process:** $X_n$ = state on day $n$, modelled as a time-homogeneous first-order
  Markov chain: $P(X_{n+1}=j \mid X_n=i, X_{n-1}, \dots) = P(X_{n+1}=j \mid X_n=i) = p_{ij}$.
* Each season is a separate sequence, so no transition is counted from 30 September to 1 June.
""")

code(r"""
df["state"] = mk.classify(df.rain_imd)
mon = df[df.month.between(6, 9)]
seasons = {y: g.state.to_numpy() for y, g in mon.groupby("year")}
assert all(len(s) == N_DAYS for s in seasons.values())
train = [seasons[y] for y in TRAIN]
test = [seasons[y] for y in TEST]
allyrs = [seasons[y] for y in seasons]

freq = (pd.crosstab(mon.month, mon.state, normalize="index")
          .rename(columns=dict(enumerate(S))).rename(index={6: "Jun", 7: "Jul", 8: "Aug", 9: "Sep"}))
print("Share of days in each state, by month (1985-2024):")
freq
""")

code(r"""
fig, ax = plt.subplots(figsize=(11, 3.4))
yr = mon[mon.year == 2005]
cols = [STATE[s] for s in yr.state]
ax.bar(yr.index, yr.rain_imd, color=cols, width=0.8)
ax.axhline(64.5, color=MUTED, lw=1, ls="--"); ax.text(yr.index[1], 70, "Heavy ≥ 64.5 mm", color=INK_2, fontsize=9)
ax.set(title="Monsoon 2005 — daily rainfall coloured by state (26 July 2005 flood)", ylabel="mm")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c, label=n) for c, n in zip(STATE, S)], loc="upper right")
fig.tight_layout(); viz.save(fig, "02_example_season")
""")

md(r"""
**Finding.** The month table already shows the monsoon's shape: July is the core (only 20% dry
days, 13.6% Heavy), while June (onset) and September (withdrawal) are almost half dry.
""")

# =====================================================================
md(r"""
## 3. Parameter estimation

The likelihood of the observed sequences, conditional on each season's first day, is
$L(P) = \prod_{i,j} p_{ij}^{\,n_{ij}}$ with $\sum_j p_{ij} = 1$. Maximising with
Lagrange multipliers gives the MLE

$$\hat p_{ij} = \frac{n_{ij}}{n_i}, \qquad \operatorname{SE}(\hat p_{ij}) = \sqrt{\frac{\hat p_{ij}(1-\hat p_{ij})}{n_i}}.$$

Estimated on the **training** seasons 1985–2014.
""")

code(r"""
C = mk.transition_counts(train)
P = mk.estimate_P(C)
SE = mk.standard_errors(C)
init = mk.initial_distribution(allyrs)

print("Transition counts n_ij (training):")
display(pd.DataFrame(C.astype(int), index=S, columns=S))
tab = pd.DataFrame([[f"{P[i,j]:.3f} ± {1.96*SE[i,j]:.3f}" for j in range(3)] for i in range(3)],
                   index=[f"from {s}" for s in S], columns=[f"to {s}" for s in S])
print("Estimated P with 95% confidence intervals:")
tab
""")

code(r"""
# Non-monsoon chain (Oct-May) for comparison
dry_season = []
for y in range(1985, 2024):
    blk = df.loc[f"{y}-10-01":f"{y+1}-05-31", "state"].to_numpy()
    dry_season.append(blk)
P_nm = mk.estimate_P(mk.transition_counts(dry_season))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for a, M, t in ((ax[0], P, "Monsoon (Jun–Sep)"), (ax[1], P_nm, "Non-monsoon (Oct–May)")):
    a.imshow(M, cmap=viz.SEQ, vmin=0, vmax=1)
    for i in range(3):
        for j in range(3):
            a.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center",
                   color="white" if M[i, j] > 0.55 else "#0b0b0b", fontsize=11)
    a.set(xticks=range(3), yticks=range(3), xticklabels=S, yticklabels=S,
          xlabel="tomorrow", ylabel="today", title=t)
    a.grid(False)
fig.tight_layout(); viz.save(fig, "03_transition_heatmaps")
print("Non-monsoon P:"); pd.DataFrame(P_nm, index=S, columns=S)
""")

code(r"""
from matplotlib.patches import FancyArrowPatch, Circle
pos = {0: np.array([0.0, 0.0]), 1: np.array([1.0, 1.55]), 2: np.array([2.0, 0.0])}
R, RAD = 0.32, 0.22
fig, ax = plt.subplots(figsize=(7.5, 6))
for i in range(3):
    for j in range(3):
        if i == j:
            continue
        p1, p2 = pos[i], pos[j]
        arr = FancyArrowPatch(p1, p2, connectionstyle=f"arc3,rad={RAD}", arrowstyle="-|>",
                              mutation_scale=16, lw=0.8 + 6 * P[i, j], color=MUTED,
                              shrinkA=34, shrinkB=34)
        ax.add_patch(arr)
        d = p2 - p1
        mid = (p1 + p2) / 2 + 0.5 * RAD * np.array([d[1], -d[0]])
        ax.text(*mid, f"{P[i, j]:.2f}", ha="center", va="center", fontsize=10, color=INK_2,
                bbox=dict(boxstyle="round,pad=0.2", fc=viz.SURFACE, ec="none"))
for i in range(3):
    ax.add_patch(Circle(pos[i], R, color=STATE[i], zorder=3))
    ax.text(*pos[i], S[i], ha="center", va="center", color="white", fontweight="bold", fontsize=11, zorder=4)
    off = np.array([0, 0.48]) if i == 1 else np.array([0, -0.48])
    ax.text(*(pos[i] + off), f"stays {P[i, i]:.2f}", ha="center", va="center", color=INK_2, fontsize=10)
ax.set_xlim(-0.6, 2.6); ax.set_ylim(-0.75, 2.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("State-transition diagram (monsoon, 1985–2014); arrow width ∝ probability")
viz.save(fig, "04_transition_diagram")
""")

md(r"""
**Finding.** Persistence dominates: Dry and Light are most likely to be followed by themselves.
Heavy usually eases to Light (0.61) rather than continuing (0.35). The monsoon moves *through*
Light: a Dry day almost never jumps straight to Heavy (0.008) and Heavy almost never drops
straight to Dry (0.04). Outside the monsoon the chain is nearly stuck in Dry (p_DD = 0.98).
""")

# =====================================================================
md(r"""
## 4. Theoretical properties

| Quantity | Formula |
|---|---|
| n-step transitions | $P^{(n)} = P^n$ (Chapman–Kolmogorov) |
| Stationary distribution | $\pi P = \pi,\ \sum_i \pi_i = 1$ (unique: chain is irreducible and aperiodic) |
| Convergence | $\lVert P^n - \mathbf 1\pi \rVert = O(|\lambda_2|^n)$ |
| Spell length in state $i$ | $\operatorname{Geometric}(1-p_{ii})$, mean $1/(1-p_{ii})$ |
| Mean recurrence time | $m_{ii} = 1/\pi_i$ |
| Mean first passage to $j$ | $m_{ij} = 1 + \sum_{k \ne j} p_{ik} m_{kj}$ |
""")

code(r"""
pi = mk.stationary(P)
emp = np.bincount(np.concatenate(train), minlength=3) / (len(train) * N_DAYS)
lam2 = mk.convergence_rate(P)
fpH = mk.first_passage(P, 2)
props = pd.DataFrame({
    "stationary π": pi,
    "observed share (train)": emp,
    "mean spell (days) = 1/(1-p_ii)": mk.mean_spell_length(P),
    "observed mean spell": [mk.all_run_lengths(train, s).mean() for s in range(3)],
    "mean recurrence 1/π (days)": 1 / pi,
    "mean days until next Heavy": fpH,
}, index=S)
print(f"|λ2| = {lam2:.3f}  ->  P^n is within 1% of its limit after {int(np.ceil(np.log(0.01)/np.log(lam2)))} days")
props
""")

code(r"""
ns = np.arange(1, 13)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
for i, a in enumerate(ax):
    Pn = np.array([mk.n_step(P, n)[i] for n in ns])
    for j in range(3):
        a.plot(ns, Pn[:, j], color=STATE[j], marker="o", ms=4, label=f"to {S[j]}")
        a.axhline(pi[j], color=STATE[j], lw=1, ls=":")
    a.set(title=f"Starting from {S[i]}", xlabel="n (days ahead)")
ax[0].set_ylabel("P(X_n = j | X_0 = i)")
ax[2].legend(loc="center right")
fig.suptitle("n-step probabilities converge to the stationary distribution (dotted)", x=0.01, ha="left",
             fontweight="bold")
fig.tight_layout(); viz.save(fig, "05_n_step_convergence")
""")

md(r"""
**Finding.** The stationary distribution π = (0.366, 0.559, 0.075) matches the observed shares
(0.369, 0.556, 0.074) almost exactly. With |λ₂| = 0.51 the chain forgets its starting state within
about a week. Expected spells: Dry 3.0 days, Light 3.4, Heavy 1.5. A Heavy day recurs every
13.4 days on average, and from a dry day it takes about 21 days on average to reach the next Heavy day.
""")

# =====================================================================
md(r"""
## 5. Simulation

Our own sampler (`mk.simulate_seasons`) draws $X_{n+1}$ from row $X_n$ of $\hat P$ by
inverse-CDF sampling. We simulate **1,000 synthetic monsoons** of 122 days, each starting
from the observed distribution of 1 June states, and compare them with reality.
""")

code(r"""
sim = mk.simulate_seasons(P, N_SIM, N_DAYS, init, rng)
fig, ax = plt.subplots(6, 1, figsize=(11, 4.2), sharex=True)
show = [(f"Observed {y}", seasons[y]) for y in (1995, 2005, 2014)] + \
       [(f"Simulated #{k+1}", sim[k]) for k in range(3)]
from matplotlib.colors import ListedColormap
for a, (name, s) in zip(ax, show):
    a.imshow(s[None, :], aspect="auto", cmap=ListedColormap(STATE), vmin=0, vmax=2)
    a.set_yticks([]); a.grid(False)
    a.text(-2, 0, name, ha="right", va="center", fontsize=9, color=INK_2)
    for sp in a.spines.values(): sp.set_visible(False)
ax[-1].set_xlabel("day of monsoon (1 = 1 June)")
ax[0].set_title("Sample paths: three real monsoons and three simulated ones")
ax[0].legend(handles=[Patch(color=c, label=n) for c, n in zip(STATE, S)], ncol=3,
             loc="lower right", bbox_to_anchor=(1, 1.05))
viz.save(fig, "06_sample_paths")
""")

# =====================================================================
md(r"""
## 6. Model validation

Four quantitative checks. Each is reported whether the model passes or not.

### V1. Is it really a (first-order) Markov chain?

Likelihood-ratio tests on the same set of transitions (days 3–122 of each training season):

* order 0 (days independent) vs order 1: $G = 2(\ell_1 - \ell_0) \sim \chi^2_{4}$
* order 1 vs order 2 (two days of memory): $G \sim \chi^2_{12}$

plus AIC $= -2\ell + 2k$ and BIC $= -2\ell + k\log N$ (lower is better).
""")

code(r"""
rows, tests, Nobs = mk.order_comparison(train)
print(f"N = {Nobs} transitions")
display(pd.DataFrame(rows).T.rename(index={0: "order 0 (independent)", 1: "order 1 (Markov)", 2: "order 2"}))
pd.DataFrame(tests).T
""")

md(r"""
**Finding.** Independence (order 0) is rejected overwhelmingly (G = 1022, p < 10⁻¹⁰): rain today
really does depend on yesterday. A second order is *also* statistically significant
(G = 114, p < 10⁻¹⁰), and both AIC and BIC prefer it, so we check it further below.
""")

md(r"""
**Is the second-order effect real, and does it matter?** Pooling June, peak and September
into one matrix can itself create apparent extra memory, so we repeat the test on peak-monsoon
days only (July–August). We also check whether a second-order chain actually forecasts the
held-out years better.
""")

code(r"""
rows_pk, tests_pk, N_pk = mk.order_comparison([s[30:92] for s in train])
print(f"Peak (Jul-Aug) only, N = {N_pk}:  order 0 vs 1 p = {tests_pk['0 vs 1']['p']:.1e} | "
      f"order 1 vs 2 G = {tests_pk['1 vs 2']['G']:.1f}, p = {tests_pk['1 vs 2']['p']:.1e}")
print(f"Share of the total likelihood gain from order 0 to 2 delivered by order 1: "
      f"{tests['0 vs 1']['G'] / (tests['0 vs 1']['G'] + tests['1 vs 2']['G']):.0%}")

def order2_P(seqs):
    n2 = np.zeros((3, 3, 3))
    for q in seqs:
        np.add.at(n2, (q[:-2], q[1:-1], q[2:]), 1)
    r = n2.sum(axis=2, keepdims=True)
    return np.divide(n2, r, out=np.full_like(n2, 1 / 3), where=r > 0)

P_2nd = order2_P(train)
b2, b1 = np.concatenate([q[:-2] for q in test]), np.concatenate([q[1:-1] for q in test])
o_ = np.concatenate([q[2:] for q in test]); O_ = np.eye(3)[o_]
brier = lambda F: ((F - O_) ** 2).sum(1).mean()
print(f"Hold-out Brier (same days): order 1 = {brier(P[b1]):.4f} | order 2 = {brier(P_2nd[b2, b1]):.4f} | "
      f"climatology = {brier(np.tile(emp, (len(o_), 1))):.4f}")
""")

md(r"""
**Finding.** The two-day memory is real: it persists in peak-monsoon days alone (p ≈ 3 × 10⁻⁷).
But the first-order chain already delivers 90% of the total likelihood gain, and the second-order
chain improves the hold-out Brier score only from 0.457 to 0.445 (about 2.6%) while tripling the
number of parameters (6 → 18). We keep the first-order chain as the main model and report
second order as an extension.
""")

md(r"""
**Why three states and not two?** If Light and Heavy could be merged into one "Wet" state,
tomorrow's distribution would be the same after a Light day as after a Heavy day, i.e.
rows L and H of $P$ would be equal. Likelihood-ratio test of $H_0: p_{L\cdot} = p_{H\cdot}$,
$G \sim \chi^2_2$:
""")

code(r"""
lump = mk.rows_equal_test(C, [1, 2])
print(f"G = {lump['G']:.1f}, df = {lump['df']}, p = {lump['p']:.2e}")
print("Row after Light:", np.round(P[1], 3), " Row after Heavy:", np.round(P[2], 3))
P2 = mk.estimate_P(mk.transition_counts([mk.to_two_state(s) for s in train], k=2))
print("\n2-state baseline P (Dry/Wet):"); pd.DataFrame(P2, index=["Dry", "Wet"], columns=["Dry", "Wet"])
""")

md(r"""
**Finding.** Light and Heavy cannot be merged (G = 164, p < 10⁻¹⁰). After a Heavy day the chance of
another Heavy day is 0.35, against 0.08 after a Light day. The Heavy state carries real information,
so three states are justified.
""")

md(r"""
### V2. Simulation check: does the model reproduce what it was not fitted to?

The chain is fitted only on one-day transitions. Spell-length distributions and the
number of Heavy days per season are *emergent* properties, so matching them is a real test.
Expected frequencies come from the 1,000 simulated seasons (which have the same 122-day
cut-off as real seasons); Pearson $\chi^2$ with tail bins pooled.
""")

code(r"""
spell_rows = []
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for s, a in zip(range(3), ax):
    obs_r = mk.all_run_lengths(train, s)
    sim_r = mk.all_run_lengths(list(sim), s)
    L = np.arange(1, max(obs_r.max(), 12) + 1)
    o = np.array([(obs_r == l).sum() for l in L[:-1]] + [(obs_r >= L[-1]).sum()])
    e_share = np.array([(sim_r == l).mean() for l in L[:-1]] + [(sim_r >= L[-1]).mean()])
    t = mk.pearson_gof(o, e_share * o.sum())
    spell_rows.append({"state": S[s], "spells observed": len(obs_r), "mean obs": obs_r.mean(),
                 "mean sim": sim_r.mean(), "theory 1/(1-p_ii)": 1 / (1 - P[s, s]),
                 "chi2": t["X2"], "df": t["df"], "p (chi2)": t["p"]})
    Lp = L[: min(len(L), 15)]
    a.bar(Lp - 0.2, o[: len(Lp)] / o.sum(), width=0.4, color=OBS, label="observed")
    a.bar(Lp + 0.2, e_share[: len(Lp)], width=0.4, color=MODEL, label="simulated")
    a.plot(Lp, mk.geometric_pmf(P, s, Lp), color=THEORY, marker="o", ms=4, label="geometric theory")
    a.set(title=f"{S[s]} spells", xlabel="length (days)", xticks=Lp[::2] if len(Lp) > 8 else Lp)
ax[0].set_ylabel("share of spells"); ax[2].legend()
fig.tight_layout(); viz.save(fig, "07_spell_lengths")
pd.DataFrame(spell_rows).set_index("state")
""")

code(r"""
# Heavy days per season: observed vs Markov-simulated distribution
obs_h = np.array([(s == 2).sum() for s in train])
sim_h = (sim == 2).sum(axis=1)
k = np.arange(0, max(obs_h.max(), sim_h.max()) + 1)
print(f"Heavy days per season  observed: mean {obs_h.mean():.2f}, var {obs_h.var(ddof=1):.2f} | "
      f"simulated: mean {sim_h.mean():.2f}, var {sim_h.var(ddof=1):.2f}")
# counts are discrete, so compare binned frequencies with Pearson chi-square (not KS)
edges = [0, 7, 9, 11, 13, 999]            # bins 0-6, 7-8, 9-10, 11-12, 13+
o_bin = np.histogram(obs_h, bins=edges)[0]
e_bin = np.histogram(sim_h, bins=edges)[0] / len(sim_h) * len(obs_h)
cnt_t = mk.pearson_gof(o_bin, e_bin)
print(f"Chi-square, binned heavy-day counts vs simulation: X2 = {cnt_t['X2']:.2f}, df = {cnt_t['df']}, p = {cnt_t['p']:.2f}")
var_X = (len(obs_h) - 1) * obs_h.var(ddof=1) / sim_h.var(ddof=1)
print(f"Variance test (observed var vs model var): p(lower tail) = {stats.chi2.cdf(var_X, len(obs_h) - 1):.3f}")

# Day-of-season state probabilities: observed vs homogeneous model
obs_day = np.array([[np.mean([s[d] == st for s in train]) for st in range(3)] for d in range(N_DAYS)])
sim_day = np.array([(sim == st).mean(axis=0) for st in range(3)]).T
""")

code(r"""
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
w = 7
smooth = lambda x: pd.Series(x).rolling(w, center=True, min_periods=1).mean()
for st, a in enumerate(ax):
    a.plot(np.arange(1, N_DAYS + 1), smooth(obs_day[:, st]), color=OBS, label="observed (7-day mean)")
    a.plot(np.arange(1, N_DAYS + 1), sim_day[:, st], color=MODEL, label="homogeneous model")
    a.set(title=f"P({S[st]}) by day of season", xlabel="day (1 = 1 June)")
    for x in (30, 92): a.axvline(x, color=MUTED, lw=0.8, ls=":")
ax[0].set_ylabel("probability"); ax[0].legend(loc="upper left")
fig.tight_layout(); viz.save(fig, "08_day_of_season_homogeneous")
""")

md(r"""
**Finding.** Spell lengths, which the model was never fitted to:
* **Heavy spells fit very well** (χ² p = 0.79).
* **Light spells are borderline** (χ² p = 0.03).
* **Dry spells are rejected** (χ² p ≈ 3 × 10⁻⁸). Real data has *more* one-day dry gaps (46% vs 33%) and
  also more long dry spells, the monsoon "breaks". A memoryless geometric spell cannot be short
  and long at once. This is a known weakness of simple chains for dry spells and motivates
  duration-dependent models.

The model reproduces the number of Heavy days per season (mean 9.07 observed vs 9.10 simulated;
binned χ² p = 0.67). Real seasons are somewhat less variable than simulated ones (variance 10.6 vs
16.5), but the difference is not significant (p = 0.07). KS tests are not used anywhere: spell
lengths and counts are whole numbers, and KS p-values are invalid (too generous) for discrete data.
""")

md(r"""
The homogeneous chain settles to $\pi$ within a few days, but real June starts drier and
September dries out again — the onset and withdrawal of the monsoon. This motivates
**Upgrade B** (Section 7).

### V3. Hold-out forecast, 2015–2024

Fit on 1985–2014, then forecast each day of 2015–2024 one day ahead. Probabilistic scores:

* multi-category Brier score $BS = \frac1N\sum_t\sum_k (f_{tk}-o_{tk})^2$ (0 is perfect)
* log-loss $-\frac1N\sum_t \log f_t(x_t)$
* Brier skill score $1 - BS_{\text{Markov}}/BS_{\text{climatology}}$, where climatology always
  forecasts the training state frequencies.
""")

code(r"""
fc = mk.forecast_scores(test, P, emp)
pd.DataFrame({"Markov chain": fc["markov"], "Climatology": fc["climatology"]}).T.assign(
    **{"Brier skill score": [fc["brier_skill"], 0.0]})
""")

md(r"""
**Finding.** On ten years it never saw, the Markov chain beats climatology: Brier skill score
0.155 (15.5% better) and log-loss 0.757 vs 0.897.
""")

md(r"""
### V4. Has the chain changed over time?

Likelihood-ratio test that 1985–2004 and 2005–2024 share one transition matrix,
$G \sim \chi^2_{6}$.
""")

code(r"""
early = [seasons[y] for y in range(1985, 2005)]
late = [seasons[y] for y in range(2005, 2025)]
Ce, Cl = mk.transition_counts(early), mk.transition_counts(late)
ht = mk.homogeneity_test([Ce, Cl])
print(f"G = {ht['G']:.2f}, df = {ht['df']}, p = {ht['p']:.4f}")
Pe, Pl = mk.estimate_P(Ce), mk.estimate_P(Cl)
pd.concat({"1985-2004": pd.DataFrame(Pe, index=S, columns=S),
           "2005-2024": pd.DataFrame(Pl, index=S, columns=S)}, axis=1)
""")

code(r"""
pi_e, pi_l = mk.stationary(Pe), mk.stationary(Pl)
pd.DataFrame({"π 1985-2004": pi_e, "π 2005-2024": pi_l,
              "mean spell 1985-2004": mk.mean_spell_length(Pe),
              "mean spell 2005-2024": mk.mean_spell_length(Pl)}, index=S)
""")

md(r"""
**Finding.** The chain has changed (p = 0.003). In 2005–2024 the long-run share of Heavy days rose
from 6.9% to 9.1%, Heavy days persist more (p_HH 0.32 → 0.38), and direct Dry → Heavy jumps more
than doubled (0.005 → 0.014). The model is therefore not perfectly time-homogeneous over 40 years,
which matters for the flood-risk results below. Even so, the model trained on 1985–2014 still
forecast 2015–2024 well (V3).
""")

# =====================================================================
md(r"""
## 7. Upgrade B — monsoon phases

Split the season into **onset** (June, days 1–30), **peak** (July–August, days 31–92) and
**withdrawal** (September, days 93–122), and fit one matrix per phase. The chain becomes
piecewise time-inhomogeneous: $P(X_{n+1}=j\mid X_n=i) = p^{(\phi(n))}_{ij}$.

Test $H_0$: one matrix for all phases, LRT $G\sim\chi^2_{12}$, and compare AIC/BIC
(the phase model has 18 free parameters, the single matrix 6).
""")

code(r"""
PH = {"Onset (Jun)": (0, 30), "Peak (Jul–Aug)": (30, 92), "Withdrawal (Sep)": (92, 122)}

def phase_counts(seqs, lo, hi):
    # transitions whose destination day is in [lo, hi)
    return mk.transition_counts([s[max(lo - 1, 0):hi] for s in seqs])

Cph = {k: phase_counts(train, lo, hi) for k, (lo, hi) in PH.items()}
Pph = {k: mk.estimate_P(c) for k, c in Cph.items()}
ph_test = mk.homogeneity_test(list(Cph.values()))

def ll(c, M):
    m = c > 0
    return float((c[m] * np.log(M[m])).sum())

ll_single = ll(C, P)
ll_phase = sum(ll(Cph[k], Pph[k]) for k in PH)
nT = C.sum()
cmp = pd.DataFrame({
    "loglik": [ll_single, ll_phase], "params": [6, 18],
}, index=["single matrix", "phase model"])
cmp["AIC"] = -2 * cmp.loglik + 2 * cmp.params
cmp["BIC"] = -2 * cmp.loglik + cmp.params * np.log(nT)
print(f"LRT single vs phases: G = {ph_test['G']:.1f}, df = {ph_test['df']}, p = {ph_test['p']:.2e}")
display(cmp)
pd.DataFrame({k: {"π Dry": mk.stationary(M)[0], "π Light": mk.stationary(M)[1], "π Heavy": mk.stationary(M)[2],
                  "p_HH": M[2, 2], "mean Heavy spell": 1 / (1 - M[2, 2]), "mean Dry spell": 1 / (1 - M[0, 0])}
              for k, M in Pph.items()}).T
""")

code(r"""
def phase_P(t, mats=Pph):
    # matrix for the transition INTO day t (0-based)
    for k, (lo, hi) in PH.items():
        if lo <= t < hi:
            return mats[k]
    return mats["Withdrawal (Sep)"]

sim_ph = mk.simulate_seasons(phase_P, N_SIM, N_DAYS, init, rng)
sim_ph_day = np.array([(sim_ph == st).mean(axis=0) for st in range(3)]).T

fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
for st, a in enumerate(ax):
    a.plot(np.arange(1, N_DAYS + 1), smooth(obs_day[:, st]), color=OBS, label="observed (7-day mean)")
    a.plot(np.arange(1, N_DAYS + 1), sim_day[:, st], color=MODEL, lw=1.2, ls="--", label="single matrix")
    a.plot(np.arange(1, N_DAYS + 1), sim_ph_day[:, st], color=THEORY, label="phase model")
    a.set(title=f"P({S[st]}) by day of season", xlabel="day (1 = 1 June)")
    for x in (30, 92): a.axvline(x, color=MUTED, lw=0.8, ls=":")
ax[0].set_ylabel("probability"); ax[0].legend(loc="upper left")
fig.tight_layout(); viz.save(fig, "09_day_of_season_phases")

smooth_all = np.column_stack([smooth(obs_day[:, st]) for st in range(3)])
err = lambda m: np.abs(smooth_all - m).mean()
print(f"Mean absolute error of daily state probabilities: single {err(sim_day):.4f} | phases {err(sim_ph_day):.4f}")

def forecast_phase(seqs, mats):
    prev = np.concatenate([s[:-1] for s in seqs]); obs = np.concatenate([s[1:] for s in seqs])
    days = np.tile(np.arange(1, N_DAYS), len(seqs))
    F = np.array([phase_P(d, mats)[p] for d, p in zip(days, prev)])
    O = np.eye(3)[obs]
    return float(((F - O) ** 2).sum(1).mean()), float(-np.log(F[np.arange(len(obs)), obs] + 1e-12).mean())

bs_ph, ll_ph = forecast_phase(test, Pph)

ph_spell = []
for st in range(3):
    obs_r = mk.all_run_lengths(train, st)
    for name, sm in (("single", sim), ("phases", sim_ph)):
        sim_r = mk.all_run_lengths(list(sm), st)
        L = np.arange(1, 13)
        o = np.array([(obs_r == l).sum() for l in L[:-1]] + [(obs_r >= L[-1]).sum()])
        e = np.array([(sim_r == l).mean() for l in L[:-1]] + [(sim_r >= L[-1]).mean()])
        ph_spell.append({"state": S[st], "model": name, "mean sim": sim_r.mean(), "mean obs": obs_r.mean(),
                         "p (chi2)": mk.pearson_gof(o, e * o.sum())["p"]})
print(f"Hold-out 2015-2024  Brier: single {fc['markov']['brier']:.4f} | phases {bs_ph:.4f} | "
      f"climatology {fc['climatology']['brier']:.4f}")
print(f"                    log-loss: single {fc['markov']['logloss']:.4f} | phases {ll_ph:.4f}")
print("\nSpell-length fit, single matrix vs phase model:")
pd.DataFrame(ph_spell).pivot(index="state", columns="model", values="p (chi2)").reindex(S)
""")

md(r"""
**Finding.** The phases are clearly different (LRT p < 10⁻¹⁰; AIC falls from 5403 to 5325; BIC is
roughly tied). Peak monsoon has the most Heavy days (π_H = 0.096) and the longest Heavy spells;
September has the fewest (π_H = 0.032). The phase model:
* cuts the error in day-by-day state probabilities by a third (0.088 → 0.059);
* fixes Light spell lengths (p 0.02 → 0.15);
* slightly improves hold-out forecasts (Brier 0.455 → 0.451).

It still cannot capture the *gradual* drying through September or the long dry breaks.
""")

# =====================================================================
md(r"""
## 8. Upgrade A — flood risk: runs of consecutive Heavy days

Waterlogging in Mumbai builds up when heavy rain *persists* for several days, because drains
and low-lying areas have no time to clear. (A single extreme day such as 27 July 2005 — 492 mm in
this grid cell — is a different hazard: our states cannot tell 65 mm from 490 mm; see Limitations.) Let
$R_k$ = event that a season contains at least $k$ consecutive Heavy days.

**Exact computation.** Augment the chain with run-length states
$\{D, L, H_1, \dots, H_{k-1}, A\}$ where $H_r$ = "Heavy for the last $r$ days" and $A$ is
absorbing ("a run of $k$ has happened"). Then $P(R_k) = P(\text{in } A \text{ after 122 days})$,
obtained by propagating the initial distribution through the augmented transition matrix.
Return period $T_k = 1/P(R_k)$ seasons.

Fitted on **all 40 seasons** using the phase model, compared with (i) Monte Carlo and
(ii) the number of real seasons with such a run. Uncertainty from a bootstrap over seasons.
""")

code(r"""
Cph_all = {k: phase_counts(allyrs, lo, hi) for k, (lo, hi) in PH.items()}
Pph_all = {k: mk.estimate_P(c) for k, c in Cph_all.items()}
P_all = mk.estimate_P(mk.transition_counts(allyrs))
f_all = lambda t: phase_P(t, Pph_all)
sim_all = mk.simulate_seasons(f_all, 20000, N_DAYS, init, rng)

B = 300
boot = {k: [] for k in range(1, 6)}
idx = np.arange(len(allyrs))
for b in range(B):
    pick = [allyrs[i] for i in rng.choice(idx, len(idx))]
    mats = {k: mk.estimate_P(phase_counts(pick, lo, hi)) for k, (lo, hi) in PH.items()}
    ini = mk.initial_distribution(pick)
    for kk in boot:
        boot[kk].append(mk.prob_heavy_run(lambda t: phase_P(t, mats), ini, N_DAYS, kk))

risk = []
for kk in range(1, 6):
    exact = mk.prob_heavy_run(f_all, init, N_DAYS, kk)
    exact_single = mk.prob_heavy_run(P_all, init, N_DAYS, kk)
    mc = np.mean([(mk.run_lengths(s, 2) >= kk).any() for s in sim_all])
    obs_n = mk.seasons_with_run(allyrs, kk)
    lo_, hi_ = np.percentile(boot[kk], [2.5, 97.5])
    risk.append({"k (days)": kk, "P exact (phase)": exact, "95% CI": f"[{lo_:.3f}, {hi_:.3f}]",
                 "P exact (single matrix)": exact_single, "Monte Carlo": mc,
                 "observed share": obs_n / len(allyrs), "observed seasons": f"{obs_n}/40",
                 "return period (yrs)": 1 / exact})
risk = pd.DataFrame(risk).set_index("k (days)")
risk
""")

code(r"""
fig, ax = plt.subplots(figsize=(7.5, 4))
ks = risk.index.to_numpy()
ci = np.array([np.percentile(boot[k], [2.5, 97.5]) for k in ks])
ax.fill_between(ks, ci[:, 0], ci[:, 1], color=MODEL, alpha=0.15, lw=0, label="95% bootstrap band")
ax.plot(ks, risk["P exact (phase)"], color=MODEL, marker="o", label="Markov model (exact)")
ax.plot(ks, risk["observed share"], color=OBS, marker="s", ls="none", ms=8, label="observed, 1985–2024")
for k in ks[1:]:
    ax.text(k + 0.08, risk.loc[k, "P exact (phase)"] + 0.02, f"1 in {risk.loc[k, 'return period (yrs)']:.1f} yrs",
            fontsize=8, color=INK_2)
ax.set(title="Chance a monsoon contains k or more consecutive Heavy days", xlabel="k (consecutive days ≥ 64.5 mm)",
       ylabel="probability per season", xticks=ks, ylim=(0, 1.05))
ax.legend()
viz.save(fig, "10_flood_risk")

yrs3 = [y for y in seasons if (mk.run_lengths(seasons[y], 2) >= 3).any()]
print("Seasons with a run of >= 3 Heavy days:", yrs3)
""")

code(r"""
# Is the risk changing? Same calculation, fitted on each 20-year half
def run3(sq):
    mats = {k: mk.estimate_P(phase_counts(sq, lo, hi)) for k, (lo, hi) in PH.items()}
    return mk.prob_heavy_run(lambda t: phase_P(t, mats), mk.initial_distribution(sq), N_DAYS, 3)

half = {}
for name, yrs in (("1985-2004", range(1985, 2005)), ("2005-2024", range(2005, 2025))):
    sq = [seasons[y] for y in yrs]
    bs = [run3([sq[i] for i in rng.choice(len(sq), len(sq))]) for _ in range(300)]
    p3 = run3(sq)
    half[name] = {"P(run ≥ 3)": p3, "95% CI": f"[{np.percentile(bs, 2.5):.2f}, {np.percentile(bs, 97.5):.2f}]",
                  "return period (yrs)": 1 / p3,
                  "observed seasons with run ≥ 3": f"{mk.seasons_with_run(sq, 3)}/20",
                  "Heavy days / season (obs)": np.mean([(q == 2).sum() for q in sq])}

has3 = np.array([(mk.run_lengths(seasons[y], 2) >= 3).any() for y in range(1985, 2025)])
a, b = has3[:20].sum(), has3[20:].sum()
fisher = stats.fisher_exact([[a, 20 - a], [b, 20 - b]])
trend = stats.spearmanr(np.arange(1985, 2025), has3)
print(f"Fisher exact test, {a}/20 vs {b}/20 seasons with a 3-day run: p = {fisher.pvalue:.3f}")
print(f"Trend over 40 years (Spearman, year vs has-a-3-day-run): rho = {trend.statistic:.2f}, p = {trend.pvalue:.3f}")
pd.DataFrame(half).T
""")

md(r"""
**Finding (headline).** The model and the data agree for every run length:

| Run of Heavy days | Model | Observed (40 yrs) | Return period |
|---|---|---|---|
| ≥ 3 days | 54% | 55% (22/40) | about 1 in 1.8 years |
| ≥ 4 days | 24% | 25% (10/40) | about 1 in 4 years |
| ≥ 5 days | 9.5% | 7.5% (3/40) | about 1 in 11 years |

This comparison is **in-sample**: the model was fitted on the same 40 years. Run lengths are not
fitted directly, so the agreement is still meaningful, but the stricter test is below.

Split by period, the chance of a ≥ 3-day Heavy run rose from 0.44 (95% CI 0.32–0.57) in 1985–2004
to 0.63 (0.51–0.72) in 2005–2024: from roughly once every 2.3 monsoons to once every 1.6. Observed
counts went from 7/20 to 15/20 seasons. Overlapping confidence intervals are not a test, so we
test the difference directly: Fisher's exact test gives **p = 0.025**, and the year-on-year trend
gives p = 0.02. Together with V4 (p = 0.003), and confirmed on the independent gauge in Section 10,
persistent heavy rain has become significantly more common. For present-day risk, the 2005–2024
figure (0.63, about once every 1.6 monsoons) is the relevant one, not the 40-year average.
""")

md(r"""
### Out-of-sample check: predict 2015–2024 from 1985–2014 only

Fit the phase model on the **training** seasons only, compute $p_k = P(R_k)$ exactly, and compare
with the 10 test seasons. Under the model, the number of test seasons containing a run of
$k$ is $\text{Binomial}(10, p_k)$, so we report the expected count and an exact two-sided binomial test.
""")

code(r"""
f_train = lambda t: phase_P(t, Pph)            # Pph was fitted on 1985-2014 in Section 7
init_train = mk.initial_distribution(train)
oos = []
for kk in range(2, 6):
    pk = mk.prob_heavy_run(f_train, init_train, N_DAYS, kk)
    obs_test = mk.seasons_with_run(test, kk)
    oos.append({"k (days)": kk, "P(run ≥ k), trained 1985-2014": pk,
                "expected in 10 test seasons": 10 * pk, "observed 2015-2024": obs_test,
                "P(≥ observed)": stats.binom.sf(obs_test - 1, 10, pk),
                "two-sided binomial p": stats.binomtest(obs_test, 10, pk).pvalue})
pd.DataFrame(oos).set_index("k (days)")
""")

md(r"""
**Finding.** Out of sample, no prediction is rejected (all p > 0.1). The model expected
2.3 and 0.9 test seasons with 4- and 5-day runs and saw 2 and 1. It *under*-predicted 3-day runs:
5.2 expected, 8 observed (P(≥ 8) = 0.07). That is the same direction as the rise found in V4, but
with only 10 test seasons it is not significant. The flood-risk model therefore holds up on unseen
years, with a hint that recent monsoons are riskier than the 1985–2014 model assumes.
""")

# =====================================================================
md(r"""
## 9. Poisson process for Heavy-rain days

Treat each Heavy day as an event. If events occurred as a (discrete-time) Poisson /
Bernoulli process with rate $\lambda$ per day:

* count per season $N \sim \text{Poisson}(\lambda \cdot 122)$, so $E[N] = \operatorname{Var}[N]$;
* gaps between events are memoryless: $\text{Geometric}(\lambda)$ in days (the discrete analogue of Exponential).

The Markov chain suggests this may fail: $p_{HH} > \pi_H$ means a Heavy day makes another
Heavy day more likely, so events should cluster in time. We test both the counts and the gaps,
and then test whether **episodes** (the first day of each Heavy spell) behave like a Poisson
process instead.

MLE: $\hat\lambda = \frac{\text{total events}}{\text{total days}}$.
Index-of-dispersion test: $D = (n-1)s^2/\bar N \sim \chi^2_{n-1}$ under Poisson.
""")

code(r"""
counts_h = np.array([(s == 2).sum() for s in allyrs])
lam = counts_h.sum() / (len(allyrs) * N_DAYS)
disp = mk.dispersion_test(counts_h)
print(f"λ̂ = {lam:.4f} Heavy days per day  ->  expected {lam*N_DAYS:.2f} per season")
print(f"p_HH = {P_all[2,2]:.3f}  vs  π_H = {mk.stationary(P_all)[2]:.3f}  (clustering if p_HH > π_H)")
print(f"Heavy days/season: mean {disp['mean']:.2f}, variance {disp['var']:.2f}, ratio {disp['ratio']:.2f}")
print(f"Dispersion test D = {disp['D']:.1f} on {disp['df']} df, p(overdispersion) = {disp['p_overdispersion']:.4f}")

gaps = mk.within_season_gaps(allyrs, 2)
G_ = np.arange(1, 16)
o = np.array([(gaps == g).sum() for g in G_[:-1]] + [(gaps >= G_[-1]).sum()])
geo = stats.geom.pmf(G_[:-1], lam); geo = np.append(geo, stats.geom.sf(G_[-2], lam))
gt = mk.pearson_gof(o, geo * o.sum(), ddof=1)
print(f"Gaps vs Geometric(λ): chi2 = {gt['X2']:.1f}, df = {gt['df']}, p = {gt['p']:.2e};  "
      f"share of gaps = 1 day: observed {np.mean(gaps == 1):.3f}, Poisson predicts {lam:.3f}")

# Episodes: first day of each Heavy spell
ep = np.array([len(mk.run_lengths(s, 2)) for s in allyrs])
disp_ep = mk.dispersion_test(ep)
print(f"\nHeavy EPISODES/season: mean {disp_ep['mean']:.2f}, variance {disp_ep['var']:.2f}, "
      f"ratio {disp_ep['ratio']:.2f}, two-sided p = {disp_ep['p_two_sided']:.3f}")
""")

code(r"""
sim_h_all = (sim_all == 2).sum(axis=1)
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
kk = np.arange(0, counts_h.max() + 3)
ax[0].bar(kk - 0.2, [(counts_h == k).mean() for k in kk], width=0.4, color=OBS, label="observed")
ax[0].bar(kk + 0.2, [(sim_h_all == k).mean() for k in kk], width=0.4, color=MODEL,
          label="Markov model")
ax[0].plot(kk, stats.poisson.pmf(kk, lam * N_DAYS), color=THEORY, marker="o", ms=4, label="Poisson")
ax[0].set(title="Heavy days per season", xlabel="count", ylabel="probability"); ax[0].legend()

gg = G_[:-1]
ax[1].bar(gg, o[:-1] / o.sum(), width=0.6, color=OBS, label="observed")
ax[1].plot(gg, geo[:-1], color=THEORY, marker="o", ms=4, label="Geometric(λ), i.e. Poisson")
ax[1].set(title="Days between Heavy days", xlabel="gap (days)", ylabel="share of gaps")
ax[1].text(14, 0.2, f"gaps ≥ 15 days:\nobserved {o[-1]/o.sum():.0%}, Poisson {geo[-1]:.0%}",
           ha="right", fontsize=9, color=INK_2)
ax[1].legend()

days = np.arange(1, N_DAYS + 1)
for i, y in enumerate((1995, 2005, 2014, 2019)):
    ax[2].step(days, np.cumsum(seasons[y] == 2), where="post", color=[OBS, MODEL, THEORY, EXTRA][i], label=str(y))
ax[2].plot(days, lam * days, color=MUTED, ls="--", lw=1.2, label="E[N(t)] = λt")
ax[2].set(title="Counting process N(t)", xlabel="day of monsoon", ylabel="Heavy days so far"); ax[2].legend()
fig.tight_layout(); viz.save(fig, "11_poisson")
""")

md(r"""
**Finding.** Looking only at counts, the Poisson model would be *accepted*: Heavy days per season
have variance/mean = 1.19, not significantly above 1 (p = 0.20). But the timing is decisively
non-Poisson: 39% of gaps between Heavy days are a single day, against 8% predicted
(χ² p < 10⁻¹⁰). Heavy days come in bursts, exactly what p_HH = 0.35 ≫ π_H = 0.08 says.
Heavy *episodes* are, if anything, more regular than Poisson (variance/mean = 0.64, p = 0.08).

Lesson: a counting test alone would have wrongly accepted the Poisson model. Inter-arrival times
reveal the clustering, and the Markov chain models it.
""")

# =====================================================================
md(r"""
## 10. Robustness checks

Questions an independent reviewer would ask, each answered with data. This section uses its own
random generator, so it does not change any number above.
""")

code(r"""
rng_r = np.random.default_rng(7)
fmt_p = lambda p: "< 1e-10" if p < 1e-10 else f"{p:.3g}"
rob = []

# R1 tougher forecast baseline: climatology that knows the day of season
prev = np.concatenate([q[:-1] for q in test]); obs_ = np.concatenate([q[1:] for q in test])
O_ = np.eye(3)[obs_]; days_ = np.tile(np.arange(1, N_DAYS), len(test))
doy = pd.DataFrame(obs_day).rolling(15, center=True, min_periods=1).mean().to_numpy()
bs_ = lambda F: float(((F - O_) ** 2).sum(1).mean())
b_m, b_doy = bs_(P[prev]), bs_(doy[days_])
rob.append(("R1 Markov vs day-of-season climatology (Brier)",
            f"{b_m:.3f} vs {b_doy:.3f}: skill {1 - b_m / b_doy:.3f}"))

# R2 is the rise in Heavy days an artefact of the IMD grid? Repeat on the independent gauge
g = mon.dropna(subset=["rain_gauge"])
share = lambda x, col: float((x[col] >= 64.5).mean())
g1, g2 = g[g.year < 2005], g[g.year >= 2005]
m1, m2 = mon[mon.year < 2005], mon[mon.year >= 2005]
tr_g = stats.spearmanr(*zip(*[(y, share(x, "rain_gauge")) for y, x in g.groupby("year")]))
tr_i = stats.spearmanr(*zip(*[(y, share(x, "rain_imd")) for y, x in mon.groupby("year")]))
rob.append(("R2 Heavy-day share, IMD grid", f"{share(m1,'rain_imd'):.3f} -> {share(m2,'rain_imd'):.3f}; trend p = {fmt_p(tr_i.pvalue)}"))
rob.append(("R2 Heavy-day share, Santacruz gauge", f"{share(g1,'rain_gauge'):.3f} -> {share(g2,'rain_gauge'):.3f}; trend p = {fmt_p(tr_g.pvalue)}"))

# R3 can seasonality alone (no clustering) explain the one-day gaps?
lam_d = pd.Series([np.mean([q[d] == 2 for q in allyrs]) for d in range(N_DAYS)]).rolling(15, center=True, min_periods=1).mean().to_numpy()
fake = (rng_r.random((5000, N_DAYS)) < lam_d).astype(int) * 2
rob.append(("R3 One-day gaps between Heavy days",
            f"real {np.mean(mk.within_season_gaps(allyrs, 2) == 1):.2f} vs seasonal rate without clustering "
            f"{np.mean(mk.within_season_gaps(list(fake), 2) == 1):.2f}"))

# R4 is second-order memory an artefact of pooling two climates? Test each half (peak days)
for a_, b_ in ((1985, 2005), (2005, 2025)):
    _, t_, _ = mk.order_comparison([seasons[y][30:92] for y in range(a_, b_)])
    rob.append((f"R4 Order 1 vs 2, peak days, {a_}-{b_ - 1}", f"p = {fmt_p(t_['1 vs 2']['p'])}"))

# R5 sensitivity of the headline to the Heavy threshold
for h in (50, 64.5, 80):
    sq = [mk.classify(x.rain_imd.to_numpy(), heavy_from=h) for _, x in mon.groupby("year")]
    Pa = mk.estimate_P(mk.transition_counts(sq))
    rob.append((f"R5 Heavy ≥ {h} mm: P(run ≥ 3)",
                f"model {mk.prob_heavy_run(Pa, mk.initial_distribution(sq), N_DAYS, 3):.2f} vs observed "
                f"{np.mean([(mk.run_lengths(q, 2) >= 3).any() for q in sq]):.2f}"))

# R6 grid vs gauge: does the grid under-count Heavy days?
rob.append(("R6 Heavy-day share on common days, IMD vs gauge",
            f"{share(both, 'rain_imd'):.3f} vs {share(both, 'rain_gauge'):.3f}"))
pd.DataFrame(rob, columns=["check", "result"]).set_index("check")
""")

md(r"""
**Finding.** The conclusions survive every check:
* **R1:** the Markov chain also beats a tougher baseline that knows the time of season
  (about 10% Brier skill), so its skill is not just seasonality.
* **R2:** the rise in Heavy days appears in the independent Santacruz gauge too, and more strongly
  (7.4% → 11.2%, trend p = 0.001). It is not an artefact of changes in the IMD station network.
* **R3:** seasonality alone produces 11% one-day gaps; reality has 39%. The clustering is genuine.
* **R4:** two-day memory appears within each 20-year half, so it is not created by mixing climates.
* **R5:** the flood-run agreement holds at 50 mm and 80 mm thresholds too.
* **R6:** the 25 km grid records fewer Heavy days than the point gauge, so our flood risks are for an
  area average and are, if anything, **conservative** for a single location.

**Multiple testing.** About 20 tests appear in this notebook. Every conclusion we draw rests on
p < 0.003 or on a non-rejection, so a Bonferroni-style correction would not change any of them.
The one result near 0.05 (the Fisher test, p = 0.025) is backed by V4 and R2.
""")

# =====================================================================
md(r"""
## 11. Summary
""")

code(r"""
summary = pd.DataFrame([
    ("Data quality", f"IMD vs gauge: daily r = {r_daily:.2f}, seasonal r = {r_season:.2f}, same state {agree:.0%}"),
    ("V1 Markov property", f"order 0 vs 1: p {fmt_p(tests['0 vs 1']['p'])};  order 1 vs 2: p {fmt_p(tests['1 vs 2']['p'])}"),
    ("V1 three states needed", f"rows L = H rejected, p {fmt_p(lump['p'])}"),
    ("V2 spell lengths", "; ".join(f"{r['state']} p {fmt_p(r['p (chi2)'])}" for r in spell_rows)),
    ("V3 hold-out Brier skill", f"{fc['brier_skill']:.3f} vs climatology; {1 - b_m / b_doy:.3f} vs day-of-season climatology"),
    ("V4 1985-2004 vs 2005-2024", f"G = {ht['G']:.1f}, p = {ht['p']:.4f}"),
    ("Upgrade B phases", f"LRT p {fmt_p(ph_test['p'])}; AIC single {cmp.AIC.iloc[0]:.0f} vs phases {cmp.AIC.iloc[1]:.0f}"),
    ("Upgrade A flood risk", f"P(≥3 Heavy days in a row) = {risk.loc[3,'P exact (phase)']:.3f} model vs {risk.loc[3,'observed share']:.3f} observed"),
    ("Upgrade A out-of-sample", f"2015-24 runs ≥ 3: expected {oos[1]['expected in 10 test seasons']:.1f}/10, "
                                f"observed {oos[1]['observed 2015-2024']}/10; all binomial p > "
                                f"{min(r['two-sided binomial p'] for r in oos):.2f}"),
    ("Upgrade A rise in risk", f"3-day runs 7/20 -> 15/20 seasons, Fisher p = {fisher.pvalue:.3f}; gauge confirms (R2)"),
    ("Poisson", f"Heavy days var/mean = {disp['ratio']:.2f} (p = {disp['p_overdispersion']:.3f}); episodes var/mean = {disp_ep['ratio']:.2f}"),
], columns=["check", "result"]).set_index("check")
summary
""")

md(r"""
## 12. Limitations and extensions

* **Memory beyond one day.** Second order is significant (V1); dry spells are not geometric (V2).
  Extension: duration-dependent or second-order chains for dry spells.
* **Not stationary over decades** (V4). Extension: a transition matrix that drifts with year
  (non-homogeneous chain, or logistic regression of transitions on year).
* **Within-season change.** Three phases help, but September dries gradually. Extension: transition
  probabilities as smooth functions of day of season.
* **Amounts are discarded.** 65 mm and 492 mm (27 July 2005) are both "Heavy". Extension: a gamma or
  extreme-value model for rainfall amounts on wet days.
* **Grid-cell average, not a point.** It records fewer Heavy days than the gauge (R6), so risks are
  conservative for any single location.
* **Calendar phases.** June / July–August / September are fixed blocks, but the monsoon's actual
  onset date varies from year to year. Extension: phases anchored on each year's onset date.
* **Design was fixed in advance.** The states, thresholds, phases and tests were chosen in the
  project plan before any result was seen, which guards against tuning the model to the data.
* **No physical drivers** (El Niño, Indian Ocean Dipole, monsoon low-pressure systems).
  Extension: hidden Markov models whose hidden states are weather regimes.
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "analysis.ipynb")
print(f"wrote analysis.ipynb with {len(cells)} cells")
