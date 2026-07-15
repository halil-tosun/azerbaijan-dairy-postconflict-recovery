"""
03_second_stage.py
===================
Simar & Wilson (2007) truncated-normal second-stage regression of
bootstrap-corrected DEA efficiency on production/economic covariates,
plus all second-stage robustness variants reported in the paper.

Reproduces: Table 4 (baseline); SI S3.1 (lagged covariates), S3.2
(bootstrap SEs), S5 (post-conflict exclusion), S6 (reduced form),
S9 (outlier-treatment sensitivity).
"""
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
import importlib.util
import os

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)
OUT_DIR = data_prep.OUT_DIR


def _negloglik_factory(y, X, upper=1.0):
    k = X.shape[1]

    def negloglik(params):
        beta = params[:k]
        sigma = np.exp(params[k])
        mu = X @ beta
        z = (y - mu) / sigma
        trunc = (upper - mu) / sigma
        return -np.sum(-np.log(sigma) + norm.logpdf(z) - norm.logcdf(trunc))
    return negloglik, k


def fit_truncreg(y, X, compute_se=True):
    negloglik, k = _negloglik_factory(y, X)
    beta0 = np.linalg.lstsq(X, y, rcond=None)[0]
    sigma0 = np.std(y - X @ beta0)
    x0 = np.concatenate([beta0, [np.log(sigma0)]])
    res1 = minimize(negloglik, x0, method='Nelder-Mead',
                     options={'maxiter': 20000, 'xatol': 1e-9, 'fatol': 1e-9})
    res = minimize(negloglik, res1.x, method='BFGS')
    beta = res.x[:k]

    se = pval = None
    if compute_se:
        eps = 1e-5
        n_par = len(res.x)
        H = np.zeros((n_par, n_par))
        for i in range(n_par):
            for j in range(n_par):
                dxi = np.zeros(n_par); dxi[i] = eps
                dxj = np.zeros(n_par); dxj[j] = eps
                H[i, j] = (negloglik(res.x + dxi + dxj) - negloglik(res.x + dxi - dxj)
                           - negloglik(res.x - dxi + dxj) + negloglik(res.x - dxi - dxj)) / (4 * eps * eps)
        cov = np.linalg.inv(H)
        se = np.sqrt(np.diag(cov))[:k]
        pval = 2 * (1 - norm.cdf(np.abs(beta / se)))
    return beta, se, pval


def build_second_stage_sample(dea_df, outlier_treatment='exclude_five'):
    d = data_prep.load_district_panel(exclude_aghdara=False, outlier_treatment=outlier_treatment)
    m = dea_df.merge(d[['region', 'year', 'cost_milk', 'profitability']], on=['region', 'year'], how='left')
    m['cattle_per1000'] = m['cows_heads'] / 1000.0
    m['fodder_per100'] = m['fodder_sown_area_ha'] / 100.0
    return m.dropna(subset=['vrs_bc', 'cattle_per1000', 'fodder_per100', 'cost_milk', 'profitability'])


def run(dea_full_csv=None):
    dea = pd.read_csv(dea_full_csv or os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'))

    # ---- Table 4: baseline ----
    d = build_second_stage_sample(dea, 'exclude_five')
    d.to_csv(os.path.join(OUT_DIR, 'table6_second_stage_sample.csv'), index=False)
    y = d['vrs_bc'].values
    X = np.column_stack([np.ones(len(d)), d['cattle_per1000'], d['fodder_per100'], d['cost_milk'], d['profitability']])
    beta, se, pval = fit_truncreg(y, X)
    names = ['Intercept', 'Cattle(1000 head)', 'Fodder(100ha)', 'Cost(AZN/centner)', 'Profitability(%)']
    table4 = pd.DataFrame({'variable': names, 'coef': beta, 'se': se, 'p_value': pval})
    table4.to_csv(os.path.join(OUT_DIR, 'table6_baseline.csv'), index=False)
    print('Table 4 (N=%d):\n' % len(d), table4)

    # ---- S3.1: lagged covariates ----
    raw = data_prep.load_district_panel(exclude_aghdara=False, outlier_treatment='exclude_five')
    raw = raw.sort_values(['region', 'year'])
    raw['cost_milk_lag'] = raw.groupby('region')['cost_milk'].shift(1)
    raw['profitability_lag'] = raw.groupby('region')['profitability'].shift(1)
    dlag = dea.merge(raw[['region', 'year', 'cost_milk_lag', 'profitability_lag']], on=['region', 'year'], how='left')
    dlag['cattle_per1000'] = dlag['cows_heads'] / 1000.0
    dlag['fodder_per100'] = dlag['fodder_sown_area_ha'] / 100.0
    dlag = dlag.dropna(subset=['vrs_bc', 'cattle_per1000', 'fodder_per100', 'cost_milk_lag', 'profitability_lag'])
    y = dlag['vrs_bc'].values
    X = np.column_stack([np.ones(len(dlag)), dlag['cattle_per1000'], dlag['fodder_per100'],
                          dlag['cost_milk_lag'], dlag['profitability_lag']])
    beta, se, pval = fit_truncreg(y, X)
    s31 = pd.DataFrame({'variable': names, 'coef': beta, 'se': se, 'p_value': pval})
    s31.to_csv(os.path.join(OUT_DIR, 's3_1_lagged.csv'), index=False)
    print('\nS3.1 lagged (N=%d):\n' % len(dlag), s31)

    # ---- S3.2: bootstrap SEs (case resampling, 300 reps) ----
    n = len(d)
    rng = np.random.default_rng(2024)
    B = 300
    y0 = d['vrs_bc'].values
    X0 = np.column_stack([np.ones(len(d)), d['cattle_per1000'], d['fodder_per100'], d['cost_milk'], d['profitability']])
    boots = np.zeros((B, X0.shape[1]))
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        bb, _, _ = fit_truncreg(y0[idx], X0[idx], compute_se=False)
        boots[b, :] = bb
    boot_se = boots.std(axis=0)
    beta0, _, _ = fit_truncreg(y0, X0, compute_se=False)
    boot_p = 2 * (1 - norm.cdf(np.abs(beta0 / boot_se)))
    s32 = pd.DataFrame({'variable': names, 'coef': beta0, 'bootstrap_se': boot_se, 'bootstrap_p': boot_p})
    s32.to_csv(os.path.join(OUT_DIR, 's3_2_bootstrap_se.csv'), index=False)
    print('\nS3.2 bootstrap SE (B=300):\n', s32)

    # ---- S5: post-conflict static comparison ----
    excl = dea[~dea['region'].isin(data_prep.POST_CONFLICT_DISTRICTS_ALL)]
    s5 = pd.DataFrame({'sample': ['All districts (baseline)', 'Excluding 10 post-conflict districts'],
                        'mean_efficiency': [dea['vrs_bc'].mean(), excl['vrs_bc'].mean()]})
    s5.to_csv(os.path.join(OUT_DIR, 's5_post_conflict_static.csv'), index=False)
    print('\nS5 static post-conflict comparison:\n', s5)

    d['post_conflict'] = d['region'].isin(data_prep.POST_CONFLICT_DISTRICTS_ALL).astype(int)
    X = np.column_stack([np.ones(len(d)), d['cattle_per1000'], d['fodder_per100'], d['post_conflict']])
    beta, se, pval = fit_truncreg(d['vrs_bc'].values, X)
    s5b = pd.DataFrame({'variable': ['Intercept', 'Cattle(1000)', 'Fodder(100ha)', 'PostConflict'],
                         'coef': beta, 'se': se, 'p_value': pval})
    s5b.to_csv(os.path.join(OUT_DIR, 's5_post_conflict_regression.csv'), index=False)

    # ---- S6: reduced form ----
    X = np.column_stack([np.ones(len(d)), d['cost_milk'], d['profitability']])
    beta, se, pval = fit_truncreg(d['vrs_bc'].values, X)
    s6 = pd.DataFrame({'variable': ['Intercept', 'Cost', 'Profitability'], 'coef': beta, 'se': se, 'p_value': pval})
    s6.to_csv(os.path.join(OUT_DIR, 's6_reduced_form.csv'), index=False)
    print('\nS6 reduced-form:\n', s6)

    # ---- S9: outlier-treatment sensitivity ----
    variants = {
        'full_sample_outliers_retained': 'retain_all',
        'excluding_gobustan_2013_only': 'exclude_gobustan_only',
        'excluding_all_five_extreme': 'exclude_five',
    }
    rows = []
    for label, treatment in variants.items():
        dv = build_second_stage_sample(dea, treatment)
        Xv = np.column_stack([np.ones(len(dv)), dv['cattle_per1000'], dv['fodder_per100'], dv['cost_milk'], dv['profitability']])
        beta, se, pval = fit_truncreg(dv['vrs_bc'].values, Xv)
        for nm, b, p in zip(names, beta, pval):
            rows.append({'variant': label, 'N': len(dv), 'variable': nm, 'coef': b, 'p_value': p})
    s9 = pd.DataFrame(rows)
    s9.to_csv(os.path.join(OUT_DIR, 's9_outlier_sensitivity.csv'), index=False)
    print('\nS9 outlier sensitivity:\n', s9)


if __name__ == '__main__':
    run()
