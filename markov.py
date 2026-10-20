"""Core functions for the Mumbai monsoon Markov chain project.

Everything here is written from scratch with numpy/scipy so the mathematics is
visible: estimation, stationary distribution, simulation, spell lengths,
first-passage times, hypothesis tests, forecast scoring and flood-run risk.

States are integers 0, 1, 2 = Dry, Light, Heavy. Rainfall data are handled as a
list of SEASONS (one integer array per year), so that no transition is ever
counted across the gap between one September and the next June.
"""

import numpy as np
from scipy import stats

STATE_NAMES = ("Dry", "Light", "Heavy")
DRY_BELOW = 2.5     # mm, IMD "rainy day" threshold
HEAVY_FROM = 64.5   # mm, IMD "heavy rainfall" threshold


# ---------------------------------------------------------------- states

def classify(rain_mm, dry_below=DRY_BELOW, heavy_from=HEAVY_FROM):
    """Map daily rainfall (mm) to states 0=Dry, 1=Light, 2=Heavy."""
    r = np.asarray(rain_mm, dtype=float)
    return np.where(r < dry_below, 0, np.where(r < heavy_from, 1, 2))


def to_two_state(seq):
    """Collapse Dry/Light/Heavy into Dry/Wet (0/1)."""
    return (np.asarray(seq) > 0).astype(int)


# ---------------------------------------------------------- estimation

def transition_counts(seasons, k=3):
    """n_ij = number of i -> j transitions, summed over all seasons."""
    n = np.zeros((k, k))
    for s in seasons:
        np.add.at(n, (s[:-1], s[1:]), 1)
    return n


def estimate_P(counts):
    """Maximum-likelihood estimate p_ij = n_ij / n_i."""
    row = counts.sum(axis=1, keepdims=True)
    return np.divide(counts, row, out=np.zeros_like(counts), where=row > 0)


def standard_errors(counts):
    """Asymptotic SE of each p_ij: sqrt(p_ij (1 - p_ij) / n_i)."""
    P = estimate_P(counts)
    n_i = counts.sum(axis=1, keepdims=True)
    return np.sqrt(P * (1 - P) / np.maximum(n_i, 1))


def initial_distribution(seasons, k=3):
    """Empirical distribution of the state on the first day of the season."""
    first = np.array([s[0] for s in seasons])
    return np.bincount(first, minlength=k) / len(first)


# ------------------------------------------------------ chain properties

def n_step(P, n):
    """n-step transition matrix P^n."""
    return np.linalg.matrix_power(P, n)


def stationary(P):
    """Solve pi P = pi, sum(pi) = 1 as a linear system."""
    k = P.shape[0]
    A = np.vstack([P.T - np.eye(k), np.ones(k)])
    b = np.zeros(k + 1)
    b[-1] = 1.0
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    return pi


def mean_recurrence(P):
    """Mean recurrence time m_ii = 1 / pi_i."""
    return 1.0 / stationary(P)


def mean_spell_length(P):
    """Spell length in state i is Geometric(1 - p_ii), mean 1 / (1 - p_ii)."""
    return 1.0 / (1.0 - np.diag(P))


def first_passage(P, target):
    """Expected days to first reach `target` from each state.

    Solves m_i = 1 + sum_{j != target} p_ij m_j, i.e. (I - Q) m = 1, where Q is P
    with the target row and column removed. Entry for the target itself is the
    mean recurrence time.
    """
    k = P.shape[0]
    others = [i for i in range(k) if i != target]
    Q = P[np.ix_(others, others)]
    m = np.linalg.solve(np.eye(len(others)) - Q, np.ones(len(others)))
    out = np.empty(k)
    out[others] = m
    out[target] = 1 + P[target, others] @ m
    return out


def convergence_rate(P):
    """Second-largest eigenvalue modulus: P^n -> 1 pi at rate |lambda_2|^n."""
    ev = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
    return ev[1]


# ------------------------------------------------------------ simulation

def simulate(P, n_days, x0, rng):
    """One sample path of length n_days starting from state x0."""
    k = P.shape[0]
    cum = np.cumsum(P, axis=1)
    x = np.empty(n_days, dtype=int)
    x[0] = x0
    u = rng.random(n_days)
    for t in range(1, n_days):
        x[t] = min(np.searchsorted(cum[x[t - 1]], u[t], side="right"), k - 1)
    return x


def simulate_seasons(P_of_day, n_seasons, n_days, init, rng):
    """Simulate many seasons.

    P_of_day is either one matrix (homogeneous chain) or a callable
    day -> matrix (phase-dependent chain). init is the day-1 distribution.
    """
    if callable(P_of_day):
        mats = [P_of_day(t) for t in range(n_days)]
    else:
        mats = [P_of_day] * n_days
    cums = [np.cumsum(M, axis=1) for M in mats]
    k = len(init)
    x0 = rng.choice(k, size=n_seasons, p=init)
    out = np.empty((n_seasons, n_days), dtype=int)
    out[:, 0] = x0
    u = rng.random((n_seasons, n_days))
    for t in range(1, n_days):
        c = cums[t][out[:, t - 1]]                     # (n_seasons, k)
        out[:, t] = np.minimum((u[:, t:t + 1] > c).sum(axis=1), k - 1)
    return out


# ---------------------------------------------------------------- spells

def run_lengths(seq, state):
    """Lengths of every maximal run of `state` in one sequence."""
    s = np.asarray(seq) == state
    if not s.any():
        return np.array([], dtype=int)
    padded = np.concatenate([[False], s, [False]]).astype(int)
    d = np.diff(padded)
    return np.flatnonzero(d == -1) - np.flatnonzero(d == 1)


def all_run_lengths(seasons, state):
    return np.concatenate([run_lengths(s, state) for s in seasons])


def geometric_pmf(P, state, lengths):
    """Theoretical spell-length pmf (1 - p) p^(L - 1), p = p_ii."""
    p = P[state, state]
    L = np.asarray(lengths)
    return (1 - p) * p ** (L - 1)


# ------------------------------------------------------- hypothesis tests

def _ll(counts, P):
    mask = counts > 0
    return float((counts[mask] * np.log(P[mask])).sum())


def order_comparison(seasons, k=3):
    """Compare order-0, order-1 and order-2 chains on the SAME transitions.

    Every model is scored on days 3..N of each season so the log-likelihoods are
    comparable; returns log-likelihood, parameter count, AIC, BIC per order and
    the two likelihood-ratio tests (0 vs 1, 1 vs 2).
    """
    n0 = np.zeros(k)
    n1 = np.zeros((k, k))
    n2 = np.zeros((k, k, k))
    for s in seasons:
        a, b, c = s[:-2], s[1:-1], s[2:]
        np.add.at(n0, c, 1)
        np.add.at(n1, (b, c), 1)
        np.add.at(n2, (a, b, c), 1)
    N = n0.sum()

    p0 = n0 / N
    ll0 = float((n0[n0 > 0] * np.log(p0[n0 > 0])).sum())
    ll1 = _ll(n1, estimate_P(n1))
    row2 = n2.sum(axis=2, keepdims=True)
    p2 = np.divide(n2, row2, out=np.zeros_like(n2), where=row2 > 0)
    m2 = n2 > 0
    ll2 = float((n2[m2] * np.log(p2[m2])).sum())

    params = {0: k - 1, 1: k * (k - 1), 2: k * k * (k - 1)}
    lls = {0: ll0, 1: ll1, 2: ll2}
    rows = {
        o: {"loglik": lls[o], "params": params[o],
            "AIC": -2 * lls[o] + 2 * params[o],
            "BIC": -2 * lls[o] + params[o] * np.log(N)}
        for o in (0, 1, 2)
    }
    tests = {}
    for lo, hi in ((0, 1), (1, 2)):
        G = 2 * (lls[hi] - lls[lo])
        df = params[hi] - params[lo]
        tests[f"{lo} vs {hi}"] = {"G": G, "df": df, "p": stats.chi2.sf(G, df)}
    return rows, tests, int(N)


def homogeneity_test(count_list):
    """LRT that several groups share one transition matrix.

    G = 2 sum_g sum_ij n_ij^g log(p_ij^g / p_ij^pooled), df = (G-1) * k(k-1)
    (minus any rows with no data).
    """
    pooled = estimate_P(sum(count_list))
    G = 0.0
    for c in count_list:
        Pg = estimate_P(c)
        m = c > 0
        G += 2 * float((c[m] * np.log(Pg[m] / pooled[m])).sum())
    k = pooled.shape[0]
    df = (len(count_list) - 1) * k * (k - 1)
    return {"G": G, "df": df, "p": stats.chi2.sf(G, df)}


def rows_equal_test(counts, rows):
    """LRT that the given rows of P are equal (tomorrow does not depend on
    which of these states today is in). Rejecting it for rows Light and Heavy
    means the two cannot be merged into one 'Wet' state.

    G = 2 sum_{i in rows} sum_j n_ij log(p_ij / p_j^pooled), df = (r-1)(k-1).
    """
    sub = counts[list(rows)]
    pooled = sub.sum(axis=0) / sub.sum()
    P = estimate_P(sub)
    m = sub > 0
    G = 2 * float((sub[m] * np.log(P[m] / np.broadcast_to(pooled, sub.shape)[m])).sum())
    df = (len(rows) - 1) * (counts.shape[1] - 1)
    return {"G": G, "df": df, "p": stats.chi2.sf(G, df)}


def pearson_gof(observed, expected, ddof=0, min_expected=5):
    """Pearson chi-square with tail bins pooled so every expected count >= 5."""
    o = list(np.asarray(observed, float))
    e = list(np.asarray(expected, float))
    while len(e) > 1 and e[-1] < min_expected:
        tail_e, tail_o = e.pop(), o.pop()
        e[-1] += tail_e
        o[-1] += tail_o
    o, e = np.array(o), np.array(e)
    e = e * o.sum() / e.sum()
    X2 = float(((o - e) ** 2 / e).sum())
    df = len(o) - 1 - ddof
    return {"X2": X2, "df": df, "p": stats.chi2.sf(X2, df), "bins": len(o)}


# ------------------------------------------------------------ forecasting

def forecast_scores(seasons, P, climatology):
    """One-day-ahead probabilistic forecasts, Markov vs climatology.

    Multi-category Brier score  BS = mean_t sum_k (f_tk - o_tk)^2
    Log loss                    LL = -mean_t log f_t(observed)
    Skill score                 SS = 1 - BS_markov / BS_climatology
    """
    k = P.shape[0]
    prev = np.concatenate([s[:-1] for s in seasons])
    obs = np.concatenate([s[1:] for s in seasons])
    O = np.eye(k)[obs]
    F_m = P[prev]
    F_c = np.tile(climatology, (len(obs), 1))
    eps = 1e-12

    def score(F):
        return {"brier": float(((F - O) ** 2).sum(axis=1).mean()),
                "logloss": float(-np.log(F[np.arange(len(obs)), obs] + eps).mean())}

    m, c = score(F_m), score(F_c)
    return {"markov": m, "climatology": c,
            "brier_skill": 1 - m["brier"] / c["brier"],
            "n": len(obs)}


# -------------------------------------------------------- flood-run risk

def prob_heavy_run(P_of_day, init, n_days, run_len, heavy=2):
    """Exact P(at least one run of >= run_len consecutive Heavy days in a season).

    Augmented chain: the non-Heavy states, 'Heavy with current run r'
    (r = 1..run_len-1) and an absorbing 'run reached' state. Push the
    distribution forward day by day; the absorbed mass is the answer.
    P_of_day may be a matrix or a callable day -> matrix.
    """
    get = P_of_day if callable(P_of_day) else (lambda t: P_of_day)
    k = len(init)
    others = [s for s in range(k) if s != heavy]
    no = len(others)
    # vector: [others..., H_1..H_{run_len-1}, absorbed]
    size = no + (run_len - 1) + 1
    v = np.zeros(size)
    for idx, s in enumerate(others):
        v[idx] = init[s]
    if run_len == 1:
        v[-1] = init[heavy]
    else:
        v[no] = init[heavy]
    for t in range(1, n_days):
        P = get(t)
        w = np.zeros(size)
        w[-1] = v[-1]
        # from non-Heavy states
        for idx, s in enumerate(others):
            for jdx, j in enumerate(others):
                w[jdx] += v[idx] * P[s, j]
            if run_len == 1:
                w[-1] += v[idx] * P[s, heavy]
            else:
                w[no] += v[idx] * P[s, heavy]
        # from Heavy with run r
        for r in range(1, run_len):
            mass = v[no + r - 1]
            for jdx, j in enumerate(others):
                w[jdx] += mass * P[heavy, j]
            if r + 1 >= run_len:
                w[-1] += mass * P[heavy, heavy]
            else:
                w[no + r] += mass * P[heavy, heavy]
        v = w
    return float(v[-1])


def seasons_with_run(seasons, run_len, heavy=2):
    """Number of seasons that contain a Heavy run of length >= run_len."""
    return int(sum((run_lengths(s, heavy) >= run_len).any() for s in seasons))


# ---------------------------------------------------------------- Poisson

def dispersion_test(counts):
    """Index-of-dispersion test: D = (n-1) s^2 / mean ~ chi2(n-1) under Poisson."""
    c = np.asarray(counts, float)
    n, m, v = len(c), c.mean(), c.var(ddof=1)
    D = (n - 1) * v / m
    p_over = stats.chi2.sf(D, n - 1)
    p_two = 2 * min(p_over, 1 - p_over)
    return {"mean": m, "var": v, "ratio": v / m, "D": D, "df": n - 1,
            "p_overdispersion": p_over, "p_two_sided": p_two}


def within_season_gaps(seasons, state=2):
    """Days between consecutive `state` days inside the same season."""
    gaps = []
    for s in seasons:
        idx = np.flatnonzero(np.asarray(s) == state)
        gaps.extend(np.diff(idx))
    return np.array(gaps, dtype=int)
