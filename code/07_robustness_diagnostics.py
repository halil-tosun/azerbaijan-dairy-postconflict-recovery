"""
07_robustness_diagnostics.py
=============================
Remaining diagnostic tables:
    SI S4  Balanced-panel comparison (districts with all 25 years)
    SI S7  Variance Inflation Factors (cattle stock, fodder area)
    SI S8  District-level efficiency ranking (highest / lowest 5)

For the S8 ranking exercise, districts with fewer than
MIN_YEARS_FOR_RANKING (20) years of usable data are excluded --
otherwise short, sparse series (in particular Aghdara district, which
has only a single usable year) can appear spuriously at the top or
bottom of the ranking. See 00_data_prep.py docstring.
"""
import pandas as pd
import numpy as np
import importlib.util
import os

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)
OUT_DIR = data_prep.OUT_DIR


def vif_manual(y, x):
    X = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ beta
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    return 1 / (1 - r2)


def run():
    dea = pd.read_csv(os.path.join(OUT_DIR, 'table3_dea_bootstrap_full.csv'))

    # ---- S4: balanced panel ----
    counts = dea.groupby('region').size()
    balanced = counts[counts == 25].index
    bal = dea[dea['region'].isin(balanced)]
    s4 = pd.DataFrame({
        'sample': ['Full sample', 'Balanced subpanel (25 years each)'],
        'N_districts': [dea['region'].nunique(), len(balanced)],
        'N_district_years': [len(dea), len(bal)],
        'mean': [dea['vrs_bc'].mean(), bal['vrs_bc'].mean()],
        'sd': [dea['vrs_bc'].std(), bal['vrs_bc'].std()],
    })
    corr = dea.groupby('year')['vrs_bc'].mean().corr(bal.groupby('year')['vrs_bc'].mean())
    print('S4 balanced panel:\n', s4)
    print('Correlation of full vs. balanced annual mean series:', corr)
    s4.to_csv(os.path.join(OUT_DIR, 's4_balanced_panel.csv'), index=False)

    # ---- S7: VIF (uses the full N=1523 production sample) ----
    cattle = dea['cows_heads'].values / 1000.0
    fodder = dea['fodder_sown_area_ha'].values / 100.0
    s7 = pd.DataFrame({'variable': ['Dairy cattle stock', 'Fodder-crop sown area'],
                        'VIF': [vif_manual(cattle, fodder), vif_manual(fodder, cattle)]})
    print('\nS7 VIF:\n', s7)
    s7.to_csv(os.path.join(OUT_DIR, 's7_vif.csv'), index=False)

    # ---- S8: district ranking (needs SFA TE) ----
    sfa = pd.read_csv(os.path.join(OUT_DIR, 'table4_sfa_district_te.csv'))[['region', 'year', 'sfa_te_translog']]
    d = dea.merge(sfa, on=['region', 'year'], how='left')
    agg = d.groupby('region').agg(dea_mean=('vrs_bc', 'mean'), sfa_mean=('sfa_te_translog', 'mean'),
                                   n_years=('year', 'count')).reset_index()
    ranked = agg[agg['n_years'] >= data_prep.MIN_YEARS_FOR_RANKING]
    top5 = ranked.sort_values('dea_mean', ascending=False).head(5)
    bottom5 = ranked.sort_values('dea_mean').head(5)
    print(f'\nS8 ranking (N>={data_prep.MIN_YEARS_FOR_RANKING} years filter, {len(ranked)}/{len(agg)} districts kept):')
    print('Top 5:\n', top5)
    print('Bottom 5:\n', bottom5)
    pd.concat([top5.assign(rank_group='highest'), bottom5.assign(rank_group='lowest')]).to_csv(
        os.path.join(OUT_DIR, 's8_district_ranking.csv'), index=False)

    # ---- Table 6 synthesis ----
    t2 = pd.read_csv(os.path.join(OUT_DIR, 'table3_summary.csv'))
    t3s = pd.read_csv(os.path.join(OUT_DIR, 'table4_s1_s2_sfa_summary.csv'))
    sfa_te = pd.read_csv(os.path.join(OUT_DIR, 'table4_sfa_district_te.csv'))['sfa_te_translog']
    t4 = pd.read_csv(os.path.join(OUT_DIR, 'table6_baseline.csv'))
    t5s = pd.read_csv(os.path.join(OUT_DIR, 'table5_summary.csv'))
    print('\nSynthesis check -- Table 2 mean:', t2.loc[t2.measure == 'DEA VRS (bootstrap-corrected)', 'mean'].values[0])
    print('Table 3 (translog) mean TE:', sfa_te.mean())
    print('Table 4 cattle coef:', t4.loc[1, 'coef'], 'fodder coef:', t4.loc[2, 'coef'])


if __name__ == '__main__':
    run()
