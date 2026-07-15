"""
05_event_study.py
==================
Event-study / difference-in-differences analysis of the effect of
post-conflict reintegration (2020 Nagorno-Karabakh ceasefire) on
district-level bootstrap-corrected DEA technical efficiency.

Reproduces: Table 6, Figure 5, SI Section S10 (robustness checks and
exogeneity discussion).

Treated districts (9; Aghdara excluded -- single observation only):
Aghdam, Fuzuli, Gubadli, Jabrayil, Kalbajar, Khojavand, Lachin, Shusha,
Zangilan. Comparison group: 58 unaffected districts. Treatment date:
2020 (post = year >= 2021).

Two specifications:
  1. Event-study (leads/lags, k=-1 as reference) -- tests parallel trends
     via a joint F-test on pre-period leads.
  2. Pooled two-way fixed-effects DiD -- single summary treatment effect.
"""
import pandas as pd
import numpy as np
from scipy.stats import t as tdist, f as fdist
import importlib.util
import os

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)
OUT_DIR = data_prep.OUT_DIR

TREATED = data_prep.EVENT_STUDY_TREATED  # 9 districts, Aghdara excluded


def _fe_ols(y, X):
    beta, _, _, _ = np.linalg.lstsq(X.values, y, rcond=None)
    resid = y - X.values @ beta
    n, k = X.shape
    sigma2 = np.sum(resid ** 2) / (n - k)
    XtX_inv = np.linalg.inv(X.values.T @ X.values)
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    return beta, se, n, k, resid


def event_study(dea, treated=TREATED, min_lead_bin=-10):
    d = dea.dropna(subset=['vrs_bc']).reset_index(drop=True).copy()
    d['treated'] = d['region'].isin(treated).astype(int)
    d['event_time'] = d['year'] - 2020
    d['et_bin'] = d['event_time'].clip(lower=min_lead_bin)

    y = d['vrs_bc'].values
    dist_dummies = pd.get_dummies(d['region'], prefix='dist', drop_first=True)
    year_dummies = pd.get_dummies(d['year'], prefix='yr', drop_first=True)
    et_bins = sorted([e for e in d['et_bin'].unique() if e != -1])
    et_df = pd.DataFrame({f'et_{e}': np.where((d['treated'] == 1) & (d['et_bin'] == e), 1, 0) for e in et_bins})
    X = pd.concat([pd.Series(1.0, index=d.index, name='const'), dist_dummies.astype(float),
                   year_dummies.astype(float), et_df.astype(float)], axis=1)

    beta, se, n, k, resid = _fe_ols(y, X)
    coefs = pd.Series(beta, index=X.columns)
    ses = pd.Series(se, index=X.columns)
    et_coefs = coefs[[c for c in X.columns if c.startswith('et_')]]
    et_coefs.index = [int(i.split('_')[1]) for i in et_coefs.index]
    et_se = ses[[c for c in X.columns if c.startswith('et_')]]
    et_se.index = [int(i.split('_')[1]) for i in et_se.index]
    et_coefs[-1] = 0.0; et_se[-1] = 0.0
    et_coefs = et_coefs.sort_index(); et_se = et_se.sort_index()

    # Joint F-test on pre-period leads (parallel-trends test)
    pre_cols = [c for c in X.columns if c.startswith('et_') and int(c.split('_')[1]) < -1]
    X_restricted = X.drop(columns=pre_cols)
    beta_r, _, _, _, resid_r = _fe_ols(y, X_restricted)
    ssr_u = np.sum(resid ** 2); ssr_r = np.sum(resid_r ** 2)
    q = len(pre_cols)
    F = ((ssr_r - ssr_u) / q) / (ssr_u / (n - k))
    p_pretrend = 1 - fdist.cdf(F, q, n - k)

    return et_coefs, et_se, F, q, n - k, p_pretrend


def pooled_did(dea, treated=TREATED, exclude_years=None, exclude_districts=None):
    d = dea.dropna(subset=['vrs_bc']).copy()
    if exclude_years:
        d = d[~d['year'].isin(exclude_years)]
    if exclude_districts:
        d = d[~d['region'].isin(exclude_districts)]
    d['treated'] = d['region'].isin(treated).astype(int)
    d['post'] = (d['year'] >= 2021).astype(int)
    d['did'] = d['treated'] * d['post']

    y = d['vrs_bc'].values
    dist_dummies = pd.get_dummies(d['region'], prefix='dist', drop_first=True)
    year_dummies = pd.get_dummies(d['year'], prefix='yr', drop_first=True)
    X = pd.concat([pd.Series(1.0, index=d.index, name='const'), dist_dummies.astype(float),
                   year_dummies.astype(float), pd.Series(d['did'].values, name='DiD', index=d.index)], axis=1)
    beta, se, n, k, resid = _fe_ols(y, X)
    coefs = pd.Series(beta, index=X.columns); ses = pd.Series(se, index=X.columns)
    t_stat = coefs['DiD'] / ses['DiD']
    pval = 2 * (1 - tdist.cdf(abs(t_stat), n - k))
    return coefs['DiD'], ses['DiD'], t_stat, pval, n


def run():
    dea = pd.read_csv(os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'))

    print('=== Table 6: Pooled DiD ===')
    b, se, t_stat, pval, n = pooled_did(dea)
    print(f'DiD = {b:.4f}, SE = {se:.4f}, t = {t_stat:.2f}, p = {pval:.4f}, N = {n}')
    table6 = pd.DataFrame({'specification': ['Pooled DiD (Treated x Post)'],
                            'estimate': [b], 'se': [se], 't_stat': [t_stat], 'p_value': [pval]})

    print('\n=== Event-study coefficients + parallel-trends F-test ===')
    et_coefs, et_se, F, df1, df2, p_pretrend = event_study(dea)
    print(et_coefs)
    print(f'\nJoint F-test (pre-period leads): F({df1},{df2}) = {F:.3f}, p = {p_pretrend:.4f}')
    table6b = pd.DataFrame({'specification': ['Joint test, pre-period leads'],
                             'estimate': [np.nan], 'se': [np.nan],
                             't_stat': [f'F({df1},{df2})={F:.2f}'], 'p_value': [p_pretrend]})
    pd.concat([table6, table6b]).to_csv(os.path.join(OUT_DIR, 'table7_event_study_summary.csv'), index=False)

    event_df = pd.DataFrame({'event_time': et_coefs.index, 'coef': et_coefs.values, 'se': et_se.values})
    event_df.to_csv(os.path.join(OUT_DIR, 'figure6_event_study.csv'), index=False)

    # ---- SI S10: robustness checks ----
    print('\n=== SI S10: Alternative treatment window (exclude 2024) ===')
    b2, se2, t2, p2, n2 = pooled_did(dea, exclude_years=[2024])
    print(f'DiD = {b2:.4f}, SE = {se2:.4f}, p = {p2:.4f}, N = {n2}')

    print('\n=== SI S10: Excluding Aghdam ===')
    treated_no_aghdam = [x for x in TREATED if x != 'Aghdam district']
    b3, se3, t3, p3, n3 = pooled_did(dea, treated=treated_no_aghdam)
    print(f'DiD = {b3:.4f}, SE = {se3:.4f}, p = {p3:.4f}, N = {n3}')

    s10 = pd.DataFrame({
        'check': ['Baseline (Table 6)', 'Excluding 2024', 'Excluding Aghdam'],
        'DiD': [b, b2, b3], 'SE': [se, se2, se3], 'p_value': [pval, p2, p3],
    })
    s10.to_csv(os.path.join(OUT_DIR, 's10_event_study_robustness.csv'), index=False)
    print('\n', s10)


if __name__ == '__main__':
    run()
