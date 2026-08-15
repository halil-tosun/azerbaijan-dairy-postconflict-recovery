"""
10_gpt_review_robustness.py
=============================
Implements four additional robustness/diagnostic analyses added to the paper
in response to an external methodological review (referred to in project
notes as "the GPT review"). These were originally run as one-off scripts
outside the replication package; this script formalizes them as part of the
reproducible pipeline so that the repository and the manuscript are fully
aligned.

Reproduces (Supporting Information numbering as used in the manuscript):
  S15 - DEA-SFA rank correlation (Spearman / Kendall), Section 2.3 (main text)
  S16 - Second-stage regression with year fixed effects, Section 4.2 (Discussion)
  S17 - Control-frontier bootstrap diagnostics by year, and raw theta>1 rates
        by treatment/period, Section 2.7.1 (Methods) and Section 4.4 (Discussion)
  S18 - (a) Control-frontier DiD with Aghdara excluded entirely from the
            sample (not just the treated group), Section 4.4 (Discussion)
        (b) Wild cluster (Rademacher) bootstrap p-values for the primary and
            secondary DiD estimates, Section 4.4 (Discussion)

Requires table3_dea_bootstrap_full.csv (01), table4_sfa_district_te.csv (02),
table7_control_frontier_scores.csv and table7_merged_scores.csv (09) to
already exist in OUT_DIR. Run after 01, 02, 03, and 09.
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, kendalltau, norm
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
m01 = _load('01_dea_bootstrap')
m03 = _load('03_second_stage')
m09 = _load('09_control_frontier_event_study')
OUT_DIR = data_prep.OUT_DIR
TREATED = data_prep.EVENT_STUDY_TREATED


# ---------------------------------------------------------------------------
# S15: DEA-SFA rank correlation
# ---------------------------------------------------------------------------
def run_s15_rank_correlation():
    print('=== S15: DEA-SFA rank correlation ===', flush=True)
    sfa = pd.read_csv(os.path.join(OUT_DIR, 'table4_sfa_district_te.csv'))
    sfa_mean = sfa.groupby('region')['sfa_te_translog'].mean().reset_index()
    sfa_mean.columns = ['region', 'sfa_mean']

    dea = pd.read_csv(os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'))
    dea_mean = dea.groupby('region')['vrs_bc'].mean().reset_index()
    dea_mean.columns = ['region', 'dea_mean']

    merged = dea_mean.merge(sfa_mean, on='region', how='inner')
    rho, p_s = spearmanr(merged['dea_mean'], merged['sfa_mean'])
    tau, p_k = kendalltau(merged['dea_mean'], merged['sfa_mean'])
    print(f'  N districts = {len(merged)}')
    print(f'  Spearman rho = {rho:.4f}, p = {p_s:.4e}')
    print(f'  Kendall tau = {tau:.4f}, p = {p_k:.4e}')

    merged.to_csv(os.path.join(OUT_DIR, 's15_dea_sfa_rank_correlation_data.csv'), index=False)
    pd.DataFrame([dict(statistic='Spearman rho', value=rho, p_value=p_s, n=len(merged)),
                  dict(statistic='Kendall tau', value=tau, p_value=p_k, n=len(merged))]
                 ).to_csv(os.path.join(OUT_DIR, 's15_dea_sfa_rank_correlation_summary.csv'), index=False)


# ---------------------------------------------------------------------------
# S16: Second-stage regression with year fixed effects
# ---------------------------------------------------------------------------
def run_s16_year_fe_second_stage():
    print('\n=== S16: Second-stage regression with year fixed effects ===', flush=True)
    dea = pd.read_csv(os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'))
    d = m03.build_second_stage_sample(dea, 'exclude_five')
    print(f'  N = {len(d)}')

    names = ['Intercept', 'Cattle', 'Fodder', 'Cost', 'Profitability']

    X0 = np.column_stack([np.ones(len(d)), d['cattle_per1000'], d['fodder_per100'],
                           d['cost_milk'], d['profitability']])
    y0 = d['vrs_bc'].values
    beta0, se0, _ = m03.fit_truncreg(y0, X0, compute_se=True)

    year_dummies = pd.get_dummies(d['year'], prefix='yr', drop_first=True).astype(float)
    X1 = np.column_stack([X0, year_dummies.values])
    beta1, se1, _ = m03.fit_truncreg(y0, X1, compute_se=True)

    rows = []
    for i, n in enumerate(names):
        t0 = beta0[i] / se0[i]
        p0 = 2 * (1 - norm.cdf(abs(t0)))
        t1 = beta1[i] / se1[i]
        p1 = 2 * (1 - norm.cdf(abs(t1)))
        rows.append(dict(variable=n, coef_baseline=beta0[i], se_baseline=se0[i], p_baseline=p0,
                          coef_year_fe=beta1[i], se_year_fe=se1[i], p_year_fe=p1))
        print(f'  {n}: baseline coef={beta0[i]:.6f} (p={p0:.4f}) | '
              f'year-FE coef={beta1[i]:.6f} (p={p1:.4f})')

    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 's16_second_stage_year_fe.csv'), index=False)

    yearly_cost = data_prep.load_district_panel(exclude_aghdara=False).groupby('year')['cost_milk'].mean()
    yearly_cost.to_csv(os.path.join(OUT_DIR, 's16_nominal_cost_price_trend.csv'))
    print(f'  Nominal cost price trend: {yearly_cost.dropna().iloc[0]:.1f} '
          f'({yearly_cost.dropna().index[0]}) -> {yearly_cost.dropna().iloc[-1]:.1f} '
          f'({yearly_cost.dropna().index[-1]})')


# ---------------------------------------------------------------------------
# S17: Control-frontier bootstrap diagnostics and raw theta>1 rates
# ---------------------------------------------------------------------------
def run_s17_control_frontier_diagnostics():
    print('\n=== S17: Control-frontier diagnostics ===', flush=True)
    d = pd.read_csv(os.path.join(OUT_DIR, 'table7_control_frontier_scores.csv'))

    by_year = d.groupby('year').agg(
        n_units=('region', 'count'),
        n_with_filtered_draws=('control_frontier_degenerate_draws', lambda x: (x > 0).sum()),
        n_all_degenerate=('control_frontier_all_draws_degenerate', 'sum'),
    ).reset_index()
    by_year.to_csv(os.path.join(OUT_DIR, 's17a_control_frontier_diagnostics_by_year.csv'), index=False)
    print(f'  Total district-years with >=1 filtered draw: '
          f'{(d["control_frontier_degenerate_draws"] > 0).sum()} / {len(d)}')
    print(f'  Total district-years with all-200 fallback: '
          f'{d["control_frontier_all_draws_degenerate"].sum()}')

    treated_post = d[(d['region'].isin(TREATED)) & (d['year'] >= 2021)]
    treated_pre = d[(d['region'].isin(TREATED)) & (d['year'] < 2020)]
    control_all = d[~d['region'].isin(TREATED)]

    rows = [
        dict(group='Control districts, all years', n=len(control_all),
             n_theta_gt1=(control_all['vrs_raw_control_frontier'] > 1).sum()),
        dict(group='Treated districts, pre-2020', n=len(treated_pre),
             n_theta_gt1=(treated_pre['vrs_raw_control_frontier'] > 1).sum()),
        dict(group='Treated districts, post-2020', n=len(treated_post),
             n_theta_gt1=(treated_post['vrs_raw_control_frontier'] > 1).sum()),
    ]
    summary = pd.DataFrame(rows)
    summary['share'] = summary['n_theta_gt1'] / summary['n']
    summary.to_csv(os.path.join(OUT_DIR, 's17b_raw_theta_gt1_by_group.csv'), index=False)
    for _, r in summary.iterrows():
        print(f'  {r["group"]}: {r["n_theta_gt1"]}/{r["n"]} ({r["share"]*100:.1f}%) raw theta > 1')


# ---------------------------------------------------------------------------
# S18a: Aghdara excluded entirely from sample
# ---------------------------------------------------------------------------
def run_s18a_aghdara_sensitivity():
    print('\n=== S18a: Control-frontier DiD, Aghdara excluded entirely ===', flush=True)
    d_full = data_prep.dea_analysis_sample(data_prep.load_district_panel(exclude_aghdara=False))
    d_no_aghdara = d_full[d_full['region'] != 'Aghdara district'].reset_index(drop=True)
    print(f'  N with Aghdara: {len(d_full)}, N without: {len(d_no_aghdara)}')

    cf_no_aghdara = m09.control_frontier_scores(d_no_aghdara, seed_base=9000)
    cf_no_aghdara.to_csv(os.path.join(OUT_DIR, 's18a_control_frontier_no_aghdara.csv'), index=False)

    res = m09.fit_did_all_se(cf_no_aghdara, 'vrs_bc_control_frontier')
    rows = []
    for k, v in res.items():
        print(f'  [{k}] DiD={v["estimate"]:.4f}, SE={v["se"]:.4f}, t={v["t"]:.2f}, '
              f'p={v["p"]:.4f}, N={v["n"]}')
        rows.append(dict(se_type=k, **v))
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 's18a_aghdara_excluded_did.csv'), index=False)


# ---------------------------------------------------------------------------
# S18b: Wild cluster (Rademacher) bootstrap
# ---------------------------------------------------------------------------
def wild_cluster_bootstrap_pvalue(y, X, cluster_ids, col_idx, B=999, seed=42):
    """Wild cluster (Rademacher) bootstrap p-value for a single coefficient,
    following Cameron, Gelbach and Miller (2008) / Cameron and Miller (2015).
    Null-imposed (restricted) DGP: residuals from the model with the target
    coefficient constrained to zero are re-weighted by +/-1 at the cluster
    level to generate the bootstrap distribution of the t-statistic."""
    rng = np.random.default_rng(seed)
    Xv = X.values if hasattr(X, 'values') else X
    n, k = Xv.shape
    groups = pd.Series(cluster_ids).values
    unique_groups = np.unique(groups)
    G = len(unique_groups)

    beta_full, _, _, _ = np.linalg.lstsq(Xv, y, rcond=None)
    resid_full = y - Xv @ beta_full
    XtX_inv = np.linalg.inv(Xv.T @ Xv)
    meat = np.zeros((k, k))
    for g in unique_groups:
        idx = groups == g
        s = Xv[idx].T @ resid_full[idx]
        meat += np.outer(s, s)
    corr = (G / (G - 1)) * ((n - 1) / (n - k))
    V = XtX_inv @ meat @ XtX_inv * corr
    se_full = np.sqrt(V[col_idx, col_idx])
    t_obs = beta_full[col_idx] / se_full

    X_restricted = np.delete(Xv, col_idx, axis=1)
    beta_r, _, _, _ = np.linalg.lstsq(X_restricted, y, rcond=None)
    resid_r = y - X_restricted @ beta_r
    fitted_r = X_restricted @ beta_r

    t_boot = np.zeros(B)
    for b in range(B):
        w = rng.choice([-1, 1], size=G)
        w_map = dict(zip(unique_groups, w))
        wt = np.array([w_map[g] for g in groups])
        y_star = fitted_r + resid_r * wt
        beta_star, _, _, _ = np.linalg.lstsq(Xv, y_star, rcond=None)
        resid_star = y_star - Xv @ beta_star
        meat_s = np.zeros((k, k))
        for g in unique_groups:
            idx = groups == g
            s = Xv[idx].T @ resid_star[idx]
            meat_s += np.outer(s, s)
        V_s = XtX_inv @ meat_s @ XtX_inv * corr
        se_s = np.sqrt(V_s[col_idx, col_idx])
        t_boot[b] = beta_star[col_idx] / se_s

    p_wild = np.mean(np.abs(t_boot) >= np.abs(t_obs))
    return dict(coef=beta_full[col_idx], se_cluster=se_full, t_obs=t_obs, p_wild=p_wild, B=B, G=G)


def run_s18b_wild_cluster_bootstrap():
    print('\n=== S18b: Wild cluster (Rademacher) bootstrap ===', flush=True)
    d = pd.read_csv(os.path.join(OUT_DIR, 'table7_merged_scores.csv'))
    rows = []
    for score_col, label in [('vrs_bc_control_frontier', 'control-frontier (primary)'),
                              ('vrs_bc', 'pooled-frontier (secondary)')]:
        dd = d.dropna(subset=[score_col]).reset_index(drop=True).copy()
        dd['treated'] = dd['region'].isin(TREATED).astype(int)
        dd['post'] = (dd['year'] >= 2021).astype(int)
        dd['did'] = dd['treated'] * dd['post']
        distd = pd.get_dummies(dd['region'], prefix='d', drop_first=True).astype(float)
        yrd = pd.get_dummies(dd['year'], prefix='y', drop_first=True).astype(float)
        X = pd.concat([pd.Series(1.0, index=dd.index, name='c'), distd, yrd,
                        pd.Series(dd['did'].values, name='DiD', index=dd.index)], axis=1)
        col_idx = list(X.columns).index('DiD')
        res = wild_cluster_bootstrap_pvalue(dd[score_col].values, X, dd['region'].values,
                                             col_idx, B=999, seed=42)
        print(f'  [{label}] coef={res["coef"]:.4f}, cluster SE={res["se_cluster"]:.4f}, '
              f't={res["t_obs"]:.3f}, wild cluster p (G={res["G"]}, B={res["B"]})={res["p_wild"]:.4f}')
        rows.append(dict(dependent_variable=label, **res))
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 's18b_wild_cluster_bootstrap.csv'), index=False)


# ---------------------------------------------------------------------------
if __name__ == '__main__':
    t0 = time.time()
    run_s15_rank_correlation()
    run_s16_year_fe_second_stage()
    run_s17_control_frontier_diagnostics()
    run_s18b_wild_cluster_bootstrap()
    run_s18a_aghdara_sensitivity()  # slowest (~9 min); run last
    print(f'\nDone in {time.time()-t0:.0f}s.')
