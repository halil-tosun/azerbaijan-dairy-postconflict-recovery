"""
04_malmquist_bootstrap.py
==========================
Region-level Bootstrap Malmquist Productivity Index (Simar & Wilson,
1999), computed on the complete 13-economic-region panel (2000-2024).

Reproduces: Table 5, Figures 2-4.

KEY METHODOLOGICAL NOTE (median vs. mean bias correction): Absheron-Khizi
sits on its own-period CRS frontier (theta=1) in most years -- a genuine
boundary/corner condition given only 13 cross-sectional units. With a
MEAN-based bootstrap bias correction this occasionally produces a
degenerate near-zero bias-corrected value for one year, contaminating the
entire 24-year cumulative product. Using the MEDIAN of the bootstrap
replications for bias correction (implemented below) resolves this.

Bootstrap: B = 300 replications.
"""
import pandas as pd
import numpy as np
from scipy.optimize import linprog
import time
import importlib.util
import os
import sys

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)
OUT_DIR = data_prep.OUT_DIR
B_REPS = 300


def theta_crs(x_eval, y_eval, Xref, Yref):
    nref, m = Xref.shape
    s = Yref.shape[1]
    c = np.zeros(nref + 1); c[0] = 1.0
    A_ub, b_ub = [], []
    for k in range(m):
        row = np.zeros(nref + 1); row[0] = -x_eval[k]; row[1:] = Xref[:, k]
        A_ub.append(row); b_ub.append(0.0)
    for r in range(s):
        row = np.zeros(nref + 1); row[1:] = -Yref[:, r]
        A_ub.append(row); b_ub.append(-y_eval[r])
    bounds = [(1e-8, None)] + [(0, None)] * nref
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=bounds, method='highs')
    return res.fun if res.success else np.nan


def theta_crs_batch(X, Y, Xref, Yref):
    return np.array([theta_crs(X[i], Y[i], Xref, Yref) for i in range(X.shape[0])])


def silverman_bw(theta):
    refl = np.concatenate([theta, 2 - theta])
    sigma = np.std(refl, ddof=1)
    iqr = np.subtract(*np.percentile(refl, [75, 25]))
    A = min(sigma, iqr / 1.349) if iqr > 0 else sigma
    return 0.9 * A * (len(refl)) ** (-1 / 5)


def pseudo_ref(X, theta_raw, rng, h):
    n = len(theta_raw)
    idx = rng.integers(0, n, size=n)
    draw = theta_raw[idx] + rng.normal(0, h, size=n)
    draw = np.where(draw > 1, 2 - draw, draw)
    draw = np.clip(draw, 1e-6, None)
    return X * (theta_raw / draw)[:, None]


def run(B=B_REPS, seed_base=5000, start_idx=0, end_idx=None):
    """Can be called repeatedly with start_idx/end_idx to process in chunks
    (each pair takes ~30-90s; 24 pairs total, ~15-20 min for full run)."""
    df = data_prep.load_region_panel().sort_values(['region', 'year'])
    regions = sorted(df['region'].unique())
    years = sorted(df['year'].unique())
    if end_idx is None:
        end_idx = len(years) - 1

    data = {}
    for yr in years:
        sub = df[df['year'] == yr].set_index('region').loc[regions]
        X = sub[['cows_heads', 'fodder_sown_area_ha']].values.astype(float)
        Y = sub[['milk_production_tons']].values.astype(float)
        data[yr] = (X, Y)

    own_theta = {yr: theta_crs_batch(*data[yr], *data[yr]) for yr in years}

    OUT_CSV = os.path.join(OUT_DIR, 'table5_malmquist_region_pairs.csv')
    done_pairs = set()
    if os.path.exists(OUT_CSV):
        prev = pd.read_csv(OUT_CSV)
        done_pairs = set(zip(prev['region'], prev['t']))

    t0 = time.time()
    all_rows = []
    for idx in range(start_idx, min(end_idx, len(years) - 1)):
        t, t1 = years[idx], years[idx + 1]
        if (regions[0], t) in done_pairs:
            continue
        Xt, Yt = data[t]; Xt1, Yt1 = data[t1]
        th_tt_raw, th_t1t1_raw = own_theta[t], own_theta[t1]
        th_tt1_raw = theta_crs_batch(Xt1, Yt1, Xt, Yt)
        th_t1t_raw = theta_crs_batch(Xt, Yt, Xt1, Yt1)

        h_t, h_t1 = silverman_bw(th_tt_raw), silverman_bw(th_t1t1_raw)
        rng = np.random.default_rng(seed_base + idx)
        boot_tt1 = np.zeros((B, len(regions))); boot_t1t = np.zeros((B, len(regions)))
        boot_tt_own = np.zeros((B, len(regions))); boot_t1t1_own = np.zeros((B, len(regions)))
        for b in range(B):
            Xt_pseudo = pseudo_ref(Xt, th_tt_raw, rng, h_t)
            Xt1_pseudo = pseudo_ref(Xt1, th_t1t1_raw, rng, h_t1)
            boot_tt1[b, :] = theta_crs_batch(Xt1, Yt1, Xt_pseudo, Yt)
            boot_t1t[b, :] = theta_crs_batch(Xt, Yt, Xt1_pseudo, Yt1)
            boot_tt_own[b, :] = theta_crs_batch(Xt, Yt, Xt_pseudo, Yt)
            boot_t1t1_own[b, :] = theta_crs_batch(Xt1, Yt1, Xt1_pseudo, Yt1)

        # NOTE: median (not mean) bias correction -- see module docstring
        th_tt1_bc = np.clip(2 * th_tt1_raw - np.median(boot_tt1, axis=0), 1e-6, None)
        th_t1t_bc = np.clip(2 * th_t1t_raw - np.median(boot_t1t, axis=0), 1e-6, None)
        th_tt_bc = np.clip(2 * th_tt_raw - np.median(boot_tt_own, axis=0), 1e-6, None)
        th_t1t1_bc = np.clip(2 * th_t1t1_raw - np.median(boot_t1t1_own, axis=0), 1e-6, None)

        for i, r in enumerate(regions):
            all_rows.append({'region': r, 't': t, 't1': t1,
                              'th_tt_bc': th_tt_bc[i], 'th_t1t1_bc': th_t1t1_bc[i],
                              'th_tt1_bc': th_tt1_bc[i], 'th_t1t_bc': th_t1t_bc[i]})
        print(f'  pair {t}-{t1} done ({time.time()-t0:.0f}s elapsed)', flush=True)
        # incremental save (survives interruption)
        batch_df = pd.DataFrame(all_rows)
        if os.path.exists(OUT_CSV):
            prev = pd.read_csv(OUT_CSV)
            batch_df = pd.concat([prev, batch_df]).drop_duplicates(subset=['region', 't'])
        batch_df.to_csv(OUT_CSV, index=False)
        all_rows = []

    res = pd.read_csv(OUT_CSV)
    if len(res) < len(regions) * (len(years) - 1):
        print(f'\nIncomplete: {len(res)}/{len(regions)*(len(years)-1)} rows. Re-run to continue.')
        return res, None

    res['TEC'] = res['th_t1t1_bc'] / res['th_tt_bc']
    res['TC'] = np.sqrt((res['th_tt1_bc'] / res['th_tt_bc']) * (res['th_t1t1_bc'] / res['th_t1t_bc']))
    res['M'] = res['TEC'] * res['TC']

    cum = res.groupby('region')[['TEC', 'TC', 'M']].prod()
    cum_pct = (cum - 1) * 100
    cum_pct['N_pairs'] = res.groupby('region').size()
    cum_pct = cum_pct.rename(columns={'TEC': 'cum_TEC_pct', 'TC': 'cum_TC_pct', 'M': 'cum_TFP_pct'})
    cum_pct = cum_pct.sort_values('cum_TFP_pct', ascending=False)
    cum_pct.to_csv(os.path.join(OUT_DIR, 'table5_summary.csv'))
    print('\nTable 5 (cumulative % change 2000-2024):\n', cum_pct)

    # Figure 2 data: cumulative TFP index time series
    pairs = res.sort_values(['region', 't'])
    rows = []
    for r in regions:
        sub = pairs[pairs['region'] == r].sort_values('t')
        cum_val = 1.0
        rows.append({'region': r, 'year': years[0], 'cum_TFP': 1.0})
        for _, row in sub.iterrows():
            cum_val *= row['M']
            rows.append({'region': r, 'year': row['t1'], 'cum_TFP': cum_val})
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'figure3_regional_tfp_timeseries.csv'), index=False)

    return res, cum_pct


if __name__ == '__main__':
    start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    end = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(start_idx=start, end_idx=end)
