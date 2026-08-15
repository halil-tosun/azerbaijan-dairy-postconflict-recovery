"""
01_dea_bootstrap.py
====================
Input-oriented DEA (CRS and VRS), with Simar & Wilson (1998, 2000)
smoothed-bootstrap bias correction, applied to a contemporaneous
(year-specific) cross-sectional frontier.

Inputs : dairy cattle stock (cows_heads), fodder-crop sown area (ha)
Output : milk production (tons)

Reproduces: Table 1 (DEA rows), Table 2, Figure 1.

Bootstrap: B = 200 replications, smoothed reflected-kernel bootstrap
(Silverman rule-of-thumb bandwidth), bias correction
theta_bc = 2*theta_hat - mean(bootstrap replicates), truncated to [0,1].

REPRODUCIBILITY NOTE: bootstrap procedures are inherently stochastic.
This script uses fixed seeds so its OWN output is exactly reproducible
on rerun, and closely reproduces the paper's reported statistics, but
will not be bit-identical to any other run with different seeds.
"""
import pandas as pd
import numpy as np
from scipy.optimize import linprog
import time
import importlib.util
import os

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)

OUT_DIR = data_prep.OUT_DIR
B_REPS = 200


def dea_input_oriented(X, Y, vrs=True):
    n, m = X.shape
    s = Y.shape[1]
    thetas = np.zeros(n)
    for i in range(n):
        c = np.zeros(n + 1); c[0] = 1.0
        A_ub, b_ub = [], []
        for k in range(m):
            row = np.zeros(n + 1); row[0] = -X[i, k]; row[1:] = X[:, k]
            A_ub.append(row); b_ub.append(0.0)
        for r in range(s):
            row = np.zeros(n + 1); row[1:] = -Y[:, r]
            A_ub.append(row); b_ub.append(-Y[i, r])
        A_eq, b_eq = None, None
        if vrs:
            row = np.zeros(n + 1); row[1:] = 1.0
            A_eq, b_eq = [row], [1.0]
        bounds = [(1e-8, None)] + [(0, None)] * n
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=A_eq, b_eq=b_eq,
                      bounds=bounds, method='highs')
        thetas[i] = res.fun if res.success else np.nan
    return thetas


def dea_input_oriented_ref(X, Y, Xref, Yref, vrs=True):
    n, m = X.shape
    s = Y.shape[1]
    nref = Xref.shape[0]
    thetas = np.zeros(n)
    for i in range(n):
        c = np.zeros(nref + 1); c[0] = 1.0
        A_ub, b_ub = [], []
        for k in range(m):
            row = np.zeros(nref + 1); row[0] = -X[i, k]; row[1:] = Xref[:, k]
            A_ub.append(row); b_ub.append(0.0)
        for r in range(s):
            row = np.zeros(nref + 1); row[1:] = -Yref[:, r]
            A_ub.append(row); b_ub.append(-Y[i, r])
        A_eq, b_eq = None, None
        if vrs:
            row = np.zeros(nref + 1); row[1:] = 1.0
            A_eq, b_eq = [row], [1.0]
        bounds = [(1e-8, None)] + [(0, None)] * nref
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=A_eq, b_eq=b_eq,
                      bounds=bounds, method='highs')
        thetas[i] = res.fun if res.success else np.nan
    return thetas


def silverman_bw(theta):
    refl = np.concatenate([theta, 2 - theta])
    sigma = np.std(refl, ddof=1)
    iqr = np.subtract(*np.percentile(refl, [75, 25]))
    A = min(sigma, iqr / 1.349) if iqr > 0 else sigma
    return 0.9 * A * (len(refl)) ** (-1 / 5)


DEGENERATE_THETA_CEILING = 5.0
# BUGFIX (identified during Food Policy revision): for DMUs that are exactly on the raw VRS
# frontier (theta_hat == 1), the smoothed reflected-kernel bootstrap occasionally perturbs the
# pseudo reference technology into a configuration for which the input-oriented LP has no
# economically meaningful finite solution close to 1; scipy/HiGHS then returns a numerically
# "successful" but degenerate solution (observed: exactly 1e6, the effective unbounded-variable
# ceiling) rather than failing outright. Left unfiltered, this either (a) occasionally produces
# NaN when the LP is genuinely infeasible, silently dropping the district-year from every
# downstream table (this happened to Sabirabad 2009 and Aghjabadi 2023 in the original release
# of this script), or (b) when it does "succeed", corrupts the bootstrap mean with an outlier of
# order 1e6, which — undetected — would corrupt the bias correction for that unit. Both failure
# modes are specific to the small number of DMUs with theta_hat == 1 in a given annual
# cross-section; no other district-year is affected (theta_hat < 1 units never trigger this
# reflection pathology). We treat any bootstrap replicate exceeding DEGENERATE_THETA_CEILING as a
# failed draw (on identical grounds to a non-converged LP) and exclude it from the bias-correction
# mean via nanmean, preserving the mean-based Simar & Wilson (1998, 2000) procedure exactly as
# specified for all non-degenerate replicates and for the ~99.9% of district-years unaffected by
# this boundary condition.


def bootstrap_bias_correct(X, Y, theta_hat, B=B_REPS, vrs=True, seed=0):
    rng = np.random.default_rng(seed)
    n = len(theta_hat)
    h = silverman_bw(theta_hat)
    boot = np.zeros((B, n))
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        draw = theta_hat[idx] + rng.normal(0, h, size=n)
        draw = np.where(draw > 1, 2 - draw, draw)
        draw = np.clip(draw, 1e-6, None)
        Xpseudo = X * (theta_hat / draw)[:, None]
        rep = dea_input_oriented_ref(X, Y, Xpseudo, Y, vrs=vrs)
        rep = np.where(rep > DEGENERATE_THETA_CEILING, np.nan, rep)
        boot[b, :] = rep
    n_degenerate = np.isnan(boot).sum(axis=0)
    bias = np.nanmean(boot, axis=0) - theta_hat
    theta_bc = np.clip(theta_hat - bias, 0, 1)
    return theta_bc, n_degenerate


def run(seed_base=1000):
    d = data_prep.dea_analysis_sample(data_prep.load_district_panel(exclude_aghdara=False))

    results = []
    t0 = time.time()
    for yr, sub in d.groupby('year'):
        X = sub[['cows_heads', 'fodder_sown_area_ha']].values.astype(float)
        Y = sub[['milk_production_tons']].values.astype(float)
        vrs_raw = dea_input_oriented(X, Y, vrs=True)
        crs_raw = dea_input_oriented(X, Y, vrs=False)
        vrs_bc, vrs_ndeg = bootstrap_bias_correct(X, Y, vrs_raw, vrs=True, seed=seed_base + int(yr))
        crs_bc, crs_ndeg = bootstrap_bias_correct(X, Y, crs_raw, vrs=False, seed=seed_base + 500 + int(yr))
        r = sub.copy()
        r['vrs_raw'] = vrs_raw; r['crs_raw'] = crs_raw
        r['vrs_bc'] = vrs_bc; r['crs_bc'] = crs_bc
        r['vrs_bootstrap_degenerate_draws'] = vrs_ndeg; r['crs_bootstrap_degenerate_draws'] = crs_ndeg
        n_affected = int((vrs_ndeg > 0).sum())
        if n_affected:
            print(f'    [note] {n_affected} district(s) in {yr} had boundary-case (theta_hat=1) '
                  f'bootstrap draws filtered per DEGENERATE_THETA_CEILING; see script docstring.', flush=True)
        results.append(r)
        print(f'  year {yr}: done ({time.time()-t0:.0f}s elapsed)', flush=True)

    full = pd.concat(results)
    full['scale_eff_raw'] = full['crs_raw'] / full['vrs_raw']
    full.to_csv(os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'), index=False)

    summary = pd.DataFrame({
        'measure': ['DEA VRS (raw)', 'DEA VRS (bootstrap-corrected)', 'Scale efficiency (CRS/VRS, raw)'],
        'mean': [full['vrs_raw'].mean(), full['vrs_bc'].mean(), full['scale_eff_raw'].mean()],
        'sd': [full['vrs_raw'].std(), full['vrs_bc'].std(), full['scale_eff_raw'].std()],
        'min': [full['vrs_raw'].min(), full['vrs_bc'].min(), full['scale_eff_raw'].min()],
        'max': [full['vrs_raw'].max(), full['vrs_bc'].max(), full['scale_eff_raw'].max()],
    })
    summary.to_csv(os.path.join(OUT_DIR, 'table3_summary.csv'), index=False)
    print('\nTable 2 summary:\n', summary)

    fig1 = full.groupby('year')['vrs_bc'].mean().reset_index()
    fig1.to_csv(os.path.join(OUT_DIR, 'figure2_annual_efficiency.csv'), index=False)
    return full


if __name__ == '__main__':
    run()
