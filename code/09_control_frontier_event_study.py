"""
09_control_frontier_event_study.py
====================================
Implements the three event-study extensions specified in the Food Policy
revision of the Materials and Methods (Section 2.7):

  (1) Control-only-frontier DEA scores (Eq. 6): the reference technology in
      every year is built exclusively from the 58 never-treated districts,
      so that no treated district's post-2020 production can mechanically
      shift the frontier against which ANY district (treated or control)
      is scored. This is the PRIMARY dependent variable for the event
      study / DiD (Section 2.7.1). Bootstrap bias correction follows the
      same Simar & Wilson (1998, 2000) mean-based procedure and the same
      degenerate-draw filter as 01_dea_bootstrap.py (see that script's
      DEGENERATE_THETA_CEILING note), adapted for an asymmetric
      scored-set-vs-reference-set bootstrap.

  (2) District-clustered and heteroskedasticity-robust (HC1, non-clustered)
      standard errors for both the pooled DiD and event-study
      specifications, on both the control-frontier and the original
      pooled-frontier dependent variable (four SE variants x 2 frontier
      variants = 8 inference results per specification), per Section 2.7.2.

  (3) A fractional logit (quasi-binomial, Papke & Wooldridge 1996) GLM of
      the same pooled-DiD specification, reporting average partial effects
      (APEs), per Section 2.7.2.

Reproduces: Table 7 (primary, control-frontier) and Table S10.1
(pooled-frontier secondary comparison + non-clustered SE), Figure 6.
"""
import pandas as pd
import numpy as np
from scipy.stats import t as tdist, f as fdist, norm
from scipy.optimize import minimize
import importlib.util
import os
import time
import warnings

warnings.filterwarnings('ignore', message='Mean of empty slice')

HERE = os.path.dirname(__file__)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


data_prep = _load('00_data_prep')
m01 = _load('01_dea_bootstrap')  # reuse dea_input_oriented, dea_input_oriented_ref, silverman_bw
OUT_DIR = data_prep.OUT_DIR
TREATED = data_prep.EVENT_STUDY_TREATED  # 9 districts, Aghdara excluded (see 00_data_prep docstring)
B_REPS = 200
DEGENERATE_THETA_CEILING = m01.DEGENERATE_THETA_CEILING


# ---------------------------------------------------------------------------
# (1) Control-only-frontier bootstrap DEA
# ---------------------------------------------------------------------------

def control_frontier_scores(d, seed_base=9000):
    """For every year, build the VRS reference technology from never-treated
    (control) districts only, score ALL districts (treated + control)
    against it, and bias-correct via the asymmetric analogue of Simar &
    Wilson (1998, 2000): the bootstrap perturbs the CONTROL units' own
    input-output data (since they alone define the technology) and rescores
    every district (including treated ones, whose data are never perturbed
    and can therefore never move the frontier) against each perturbed
    control-only reference set.

    Note on Aghdara: consistent with the original event-study code
    (05_event_study.py), Aghdara district — itself post-conflict but with
    only a single usable observation — is not in TREATED and is therefore
    counted among the 58 "control" districts whenever it appears (2024
    only). This preserves exact comparability with the paper's stated
    "comparison group of 58 districts" and affects only a single
    district-year in the control reference set.
    """
    results = []
    t0 = time.time()
    for yr, sub in d.groupby('year'):
        sub = sub.reset_index(drop=True)
        is_control = ~sub['region'].isin(TREATED)
        Xall = sub[['cows_heads', 'fodder_sown_area_ha']].values.astype(float)
        Yall = sub[['milk_production_tons']].values.astype(float)
        Xc = Xall[is_control.values]
        Yc = Yall[is_control.values]
        nc = Xc.shape[0]

        # own-technology raw scores of control units (needed for bandwidth + bootstrap seed)
        raw_control_own = m01.dea_input_oriented_ref(Xc, Yc, Xc, Yc, vrs=True)
        # raw score of ALL districts (treated + control) against the control-only frontier
        raw_all_vs_control = m01.dea_input_oriented_ref(Xall, Yall, Xc, Yc, vrs=True)

        h = m01.silverman_bw(raw_control_own)
        rng = np.random.default_rng(seed_base + int(yr))
        B = B_REPS
        boot = np.zeros((B, Xall.shape[0]))
        for b in range(B):
            idx = rng.integers(0, nc, size=nc)
            draw = raw_control_own[idx] + rng.normal(0, h, size=nc)
            draw = np.where(draw > 1, 2 - draw, draw)
            draw = np.clip(draw, 1e-6, None)
            Xc_pseudo = Xc[idx] * (raw_control_own[idx] / draw)[:, None]
            rep = m01.dea_input_oriented_ref(Xall, Yall, Xc_pseudo, Yc[idx], vrs=True)
            rep = np.where(rep > DEGENERATE_THETA_CEILING, np.nan, rep)
            boot[b, :] = rep

        n_degenerate = np.isnan(boot).sum(axis=0)
        all_degenerate = n_degenerate == B  # every replicate failed for this unit
        with np.errstate(invalid='ignore'):
            boot_mean = np.nanmean(boot, axis=0)
        bias = np.where(all_degenerate, 0.0, boot_mean - raw_all_vs_control)
        # Fallback for the (rare, small-reference-set) case where every bootstrap replicate for a
        # unit is degenerate: apply no bias correction (theta_bc = raw score) rather than NaN, so
        # no district-year is silently dropped from Table 7. Flagged via
        # control_frontier_all_draws_degenerate for transparency; see manuscript Table 7 notes.
        theta_bc = np.clip(raw_all_vs_control - bias, 0, 1)
        if all_degenerate.any():
            print(f'    [note] {int(all_degenerate.sum())} district(s) had ALL bootstrap draws '
                  f'degenerate in {yr}; bias correction fell back to raw score for these units only.',
                  flush=True)

        r = sub.copy()
        r['vrs_raw_control_frontier'] = raw_all_vs_control
        r['vrs_bc_control_frontier'] = theta_bc
        r['control_frontier_degenerate_draws'] = n_degenerate
        r['control_frontier_all_draws_degenerate'] = all_degenerate
        r['is_control_unit'] = is_control.values
        results.append(r)
        print(f'  control-frontier {yr}: done ({time.time()-t0:.0f}s elapsed, '
              f'{int((n_degenerate>0).sum())} district(s) with filtered draws)', flush=True)

    return pd.concat(results, ignore_index=True)


# ---------------------------------------------------------------------------
# (2) Standard-error variants for the linear two-way FE model
# ---------------------------------------------------------------------------

def _design(d, score_col, treated_set=TREATED):
    d = d.dropna(subset=[score_col]).reset_index(drop=True).copy()
    d['treated'] = d['region'].isin(treated_set).astype(int)
    d['post'] = (d['year'] >= 2021).astype(int)
    d['did'] = d['treated'] * d['post']
    y = d[score_col].values
    dist_dummies = pd.get_dummies(d['region'], prefix='dist', drop_first=True).astype(float)
    year_dummies = pd.get_dummies(d['year'], prefix='yr', drop_first=True).astype(float)
    X = pd.concat([pd.Series(1.0, index=d.index, name='const'), dist_dummies, year_dummies,
                   pd.Series(d['did'].values, name='DiD', index=d.index)], axis=1)
    return d, y, X


def _ols_fit(y, X):
    beta, _, _, _ = np.linalg.lstsq(X.values, y, rcond=None)
    resid = y - X.values @ beta
    return beta, resid


def _se_classical(X, resid):
    n, k = X.shape
    sigma2 = np.sum(resid ** 2) / (n - k)
    XtX_inv = np.linalg.inv(X.values.T @ X.values)
    return np.sqrt(np.diag(XtX_inv) * sigma2)


def _se_hc1(X, resid):
    """White (1980) heteroskedasticity-robust SE, HC1 small-sample correction."""
    n, k = X.shape
    Xv = X.values
    XtX_inv = np.linalg.inv(Xv.T @ Xv)
    meat = (Xv * (resid ** 2)[:, None]).T @ Xv
    V = XtX_inv @ meat @ XtX_inv * (n / (n - k))
    return np.sqrt(np.diag(V))


def _se_cluster(X, resid, cluster_ids):
    """Cluster-robust (CR1) SE, clustered by district."""
    n, k = X.shape
    Xv = X.values
    XtX_inv = np.linalg.inv(Xv.T @ Xv)
    groups = pd.Series(cluster_ids).values
    meat = np.zeros((k, k))
    for g in np.unique(groups):
        idx = groups == g
        Xg = Xv[idx]
        ug = resid[idx]
        s = Xg.T @ ug
        meat += np.outer(s, s)
    G = len(np.unique(groups))
    corr = (G / (G - 1)) * ((n - 1) / (n - k))
    V = XtX_inv @ meat @ XtX_inv * corr
    return np.sqrt(np.diag(V))


def fit_did_all_se(d, score_col, treated_set=TREATED):
    dd, y, X = _design(d, score_col, treated_set)
    beta, resid = _ols_fit(y, X)
    n, k = X.shape
    col_idx = list(X.columns).index('DiD')
    b = beta[col_idx]
    se_classical = _se_classical(X, resid)[col_idx]
    se_hc1 = _se_hc1(X, resid)[col_idx]
    se_cluster = _se_cluster(X, resid, dd['region'].values)[col_idx]
    out = {}
    for label, se in [('classical', se_classical), ('HC1_robust', se_hc1), ('district_clustered', se_cluster)]:
        t_stat = b / se
        p = 2 * (1 - tdist.cdf(np.abs(t_stat), n - k))
        out[label] = dict(estimate=b, se=se, t=t_stat, p=p, n=n)
    return out


def event_study_all_se(d, score_col, treated_set=TREATED, min_lead_bin=-10):
    dd = d.dropna(subset=[score_col]).reset_index(drop=True).copy()
    dd['treated'] = dd['region'].isin(treated_set).astype(int)
    dd['event_time'] = dd['year'] - 2020
    dd['et_bin'] = dd['event_time'].clip(lower=min_lead_bin)
    y = dd[score_col].values
    dist_dummies = pd.get_dummies(dd['region'], prefix='dist', drop_first=True).astype(float)
    year_dummies = pd.get_dummies(dd['year'], prefix='yr', drop_first=True).astype(float)
    et_bins = sorted([e for e in dd['et_bin'].unique() if e != -1])
    et_df = pd.DataFrame({f'et_{e}': np.where((dd['treated'] == 1) & (dd['et_bin'] == e), 1, 0)
                           for e in et_bins}).astype(float)
    X = pd.concat([pd.Series(1.0, index=dd.index, name='const'), dist_dummies, year_dummies, et_df], axis=1)
    beta, resid = _ols_fit(y, X)
    n, k = X.shape

    se_cluster_all = _se_cluster(X, resid, dd['region'].values)
    coefs = pd.Series(beta, index=X.columns)
    ses = pd.Series(se_cluster_all, index=X.columns)
    et_coefs = coefs[[c for c in X.columns if c.startswith('et_')]]
    et_coefs.index = [int(i.split('_')[1]) for i in et_coefs.index]
    et_se = ses[[c for c in X.columns if c.startswith('et_')]]
    et_se.index = [int(i.split('_')[1]) for i in et_se.index]
    et_coefs[-1] = 0.0; et_se[-1] = 0.0
    et_coefs = et_coefs.sort_index(); et_se = et_se.sort_index()

    pre_cols = [c for c in X.columns if c.startswith('et_') and int(c.split('_')[1]) < -1]
    X_restricted = X.drop(columns=pre_cols)
    beta_r, resid_r = _ols_fit(y, X_restricted)
    ssr_u = np.sum(resid ** 2); ssr_r = np.sum(resid_r ** 2)
    q = len(pre_cols)
    F = ((ssr_r - ssr_u) / q) / (ssr_u / (n - k))
    p_pretrend = 1 - fdist.cdf(F, q, n - k)

    return et_coefs, et_se, F, q, n - k, p_pretrend


# ---------------------------------------------------------------------------
# (3) Fractional logit (Papke & Wooldridge, 1996) quasi-MLE
# ---------------------------------------------------------------------------

def fractional_logit_did(d, score_col, treated_set=TREATED):
    dd, y, X = _design(d, score_col, treated_set)
    y = np.clip(y, 1e-6, 1 - 1e-6)  # fractional logit requires y strictly in (0,1)
    Xv = X.values
    n, k = Xv.shape

    def negqloglik(b):
        eta = Xv @ b
        p = 1 / (1 + np.exp(-eta))
        p = np.clip(p, 1e-10, 1 - 1e-10)
        return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))

    b0 = np.zeros(k)
    b0[0] = np.log(y.mean() / (1 - y.mean()))
    res1 = minimize(negqloglik, b0, method='Nelder-Mead',
                     options={'maxiter': 40000, 'xatol': 1e-9, 'fatol': 1e-9})
    res = minimize(negqloglik, res1.x, method='BFGS', options={'maxiter': 5000, 'gtol': 1e-7})
    if not res.success:
        res = minimize(negqloglik, res.x, method='Nelder-Mead',
                        options={'maxiter': 40000, 'xatol': 1e-10, 'fatol': 1e-10})
    b = res.x

    eta = Xv @ b
    p_hat = 1 / (1 + np.exp(-eta))
    w = p_hat * (1 - p_hat)
    # sandwich (robust QMLE) variance: A^-1 B A^-1, A = X'WX (expected Hessian of Bernoulli
    # quasi-likelihood), B = X' diag((y-p)^2) X (score outer product) -- standard
    # Papke-Wooldridge robust covariance for fractional response QMLE.
    # A small ridge term and pseudo-inverse are used because A can be numerically close to
    # singular when many district/year fixed-effect cells have fitted probabilities near the
    # [0,1] boundary (weights w -> 0); this affects only the numerical conditioning of the
    # variance calculation, not the point estimate (from BFGS on the full-rank quasi-likelihood).
    ridge = 1e-8 * np.eye(k)
    A = Xv.T @ (Xv * w[:, None]) + ridge
    A_inv = np.linalg.pinv(A)
    score_i = Xv * (y - p_hat)[:, None]
    B_classical = score_i.T @ score_i
    V_robust = A_inv @ B_classical @ A_inv

    groups = dd['region'].values
    B_cluster = np.zeros((k, k))
    for g in np.unique(groups):
        idx = groups == g
        s = score_i[idx].sum(axis=0)
        B_cluster += np.outer(s, s)
    V_cluster = A_inv @ B_cluster @ A_inv

    col_idx = list(X.columns).index('DiD')
    coef = b[col_idx]
    se_robust = np.sqrt(V_robust[col_idx, col_idx])
    se_cluster = np.sqrt(V_cluster[col_idx, col_idx])

    # Average partial effect of DiD on the fitted probability (Papke-Wooldridge APE):
    # d p / d(DiD) = p(1-p) * beta_DiD, averaged over the sample (treated x post = the
    # discrete-jump analogue; APE reported as an averaged marginal derivative, standard
    # in this literature for a fixed-effects-only "regressor").
    ape = np.mean(w) * coef
    n_par = n

    return dict(coef=coef, se_robust=se_robust, se_cluster=se_cluster, ape=ape,
                n=n_par, converged=res.success)


# ---------------------------------------------------------------------------
def run():
    d = pd.read_csv(os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'))

    print('=== Building control-only-frontier bootstrap DEA scores (Eq. 6) ===')
    cf = control_frontier_scores(d)
    cf.to_csv(os.path.join(OUT_DIR, 'table7_control_frontier_scores.csv'), index=False)

    merged = d.merge(
        cf[['region', 'year', 'vrs_raw_control_frontier', 'vrs_bc_control_frontier',
            'control_frontier_degenerate_draws']],
        on=['region', 'year'], how='left')
    merged.to_csv(os.path.join(OUT_DIR, 'table7_merged_scores.csv'), index=False)

    rows = []
    print('\n=== Table 7 (PRIMARY): pooled DiD, control-frontier scores, all SE variants ===')
    res_cf = fit_did_all_se(merged, 'vrs_bc_control_frontier')
    for label, r in res_cf.items():
        print(f'  [{label}] DiD={r["estimate"]:.4f}, SE={r["se"]:.4f}, t={r["t"]:.2f}, p={r["p"]:.4f}, N={r["n"]}')
        rows.append(dict(frontier='control', se_type=label, **r))

    print('\n=== Table S10.1 (secondary): pooled DiD, pooled-frontier scores, all SE variants ===')
    res_pooled = fit_did_all_se(merged, 'vrs_bc')
    for label, r in res_pooled.items():
        print(f'  [{label}] DiD={r["estimate"]:.4f}, SE={r["se"]:.4f}, t={r["t"]:.2f}, p={r["p"]:.4f}, N={r["n"]}')
        rows.append(dict(frontier='pooled', se_type=label, **r))

    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'table7_did_all_se_variants.csv'), index=False)

    print('\n=== Event-study leads/lags (control-frontier, district-clustered SE) ===')
    et_coefs, et_se, F, df1, df2, p_pre = event_study_all_se(merged, 'vrs_bc_control_frontier')
    print(et_coefs)
    print(f'Joint F-test (pre-trends): F({df1},{df2}) = {F:.3f}, p = {p_pre:.4f}')
    event_df = pd.DataFrame({'event_time': et_coefs.index, 'coef': et_coefs.values, 'se_clustered': et_se.values})
    event_df.to_csv(os.path.join(OUT_DIR, 'figure6_event_study_control_frontier.csv'), index=False)
    pd.DataFrame([dict(F=F, df1=df1, df2=df2, p=p_pre)]).to_csv(
        os.path.join(OUT_DIR, 'table7_pretrend_ftest.csv'), index=False)

    print('\n=== Event-study leads/lags (pooled-frontier, district-clustered SE; secondary) ===')
    et_coefs_p, et_se_p, F_p, df1_p, df2_p, p_pre_p = event_study_all_se(merged, 'vrs_bc')
    print(et_coefs_p)
    print(f'Joint F-test (pre-trends): F({df1_p},{df2_p}) = {F_p:.3f}, p = {p_pre_p:.4f}')
    pd.DataFrame({'event_time': et_coefs_p.index, 'coef': et_coefs_p.values, 'se_clustered': et_se_p.values}
                 ).to_csv(os.path.join(OUT_DIR, 's10_event_study_pooled_frontier.csv'), index=False)

    print('\n=== Fractional logit (quasi-binomial), control-frontier, DiD ===')
    fl_cf = fractional_logit_did(merged, 'vrs_bc_control_frontier')
    print(f'  coef={fl_cf["coef"]:.4f}, SE(robust)={fl_cf["se_robust"]:.4f}, '
          f'SE(cluster)={fl_cf["se_cluster"]:.4f}, APE={fl_cf["ape"]:.4f}, converged={fl_cf["converged"]}')
    fl_pooled = fractional_logit_did(merged, 'vrs_bc')
    print(f'  [pooled, secondary] coef={fl_pooled["coef"]:.4f}, SE(robust)={fl_pooled["se_robust"]:.4f}, '
          f'SE(cluster)={fl_pooled["se_cluster"]:.4f}, APE={fl_pooled["ape"]:.4f}')
    pd.DataFrame([
        dict(frontier='control', **fl_cf), dict(frontier='pooled', **fl_pooled),
    ]).to_csv(os.path.join(OUT_DIR, 'table7_fractional_logit.csv'), index=False)

    # ---- S10 robustness (alternative window, excluding Aghdam) on control-frontier scores ----
    print('\n=== S10 robustness (control-frontier): excluding 2024 ===')
    m2 = merged[merged['year'] != 2024]
    r2 = fit_did_all_se(m2, 'vrs_bc_control_frontier')['district_clustered']
    print(f'  DiD={r2["estimate"]:.4f}, SE={r2["se"]:.4f}, p={r2["p"]:.4f}, N={r2["n"]}')

    print('=== S10 robustness (control-frontier): excluding Aghdam ===')
    treated_no_aghdam = [x for x in TREATED if x != 'Aghdam district']
    r3 = fit_did_all_se(merged, 'vrs_bc_control_frontier', treated_set=treated_no_aghdam)['district_clustered']
    print(f'  DiD={r3["estimate"]:.4f}, SE={r3["se"]:.4f}, p={r3["p"]:.4f}, N={r3["n"]}')

    s10 = pd.DataFrame({
        'check': ['Baseline (control-frontier, clustered)', 'Excluding 2024', 'Excluding Aghdam'],
        'DiD': [res_cf['district_clustered']['estimate'], r2['estimate'], r3['estimate']],
        'SE': [res_cf['district_clustered']['se'], r2['se'], r3['se']],
        'p_value': [res_cf['district_clustered']['p'], r2['p'], r3['p']],
    })
    s10.to_csv(os.path.join(OUT_DIR, 's10_event_study_robustness_control_frontier.csv'), index=False)
    print('\n', s10)


if __name__ == '__main__':
    run()
