"""
06_labor_capital_land_robustness.py
=====================================
Three independently-sourced additional-input robustness checks:
  - Labor (total labor-hours in milk production; enterprise subsector)
  - Capital (agricultural capital stock; private-farm subsector)
  - Land (total land in ownership/use; enterprise subsector)

Reproduces: Table 7 (labor), SI Table S12 (capital), SI Table S13 (land).

IMPORTANT DATA-SCOPE CAVEAT (see paper Section 3.7 / SI S12): capital and
land are each available for only ONE farm-organization subsector while
output and cattle stock are reported for all categories combined. These
two variables should be read as subsector-specific proxies, not exact
all-category totals.

A three-input stochastic frontier specification is poorly behaved (see
sfa_3input_diagnostic() below, Waldman 1982 positive-skewness problem),
so all three checks are conducted primarily via DEA.
"""
import pandas as pd
import numpy as np
from scipy.optimize import linprog, minimize
from scipy.stats import norm, chi2, pearsonr, skew
import importlib.util
import os

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)
OUT_DIR = data_prep.OUT_DIR


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


def sfa_loglik(y, X):
    n, k = X.shape

    def negloglik(params):
        beta = params[:k]; sigma2 = np.exp(params[k]); lam = params[k + 1]
        sigma = np.sqrt(sigma2); eps = y - X @ beta; z = eps / sigma
        return -np.sum(np.log(2) - np.log(sigma) + norm.logpdf(z) + norm.logcdf(-z * lam))

    beta0 = np.linalg.lstsq(X, y, rcond=None)[0]
    sigma2_0 = np.var(y - X @ beta0)
    best = None
    for lam0 in (0.5, 1, 1.5, 2, 3):
        x0 = np.concatenate([beta0, [np.log(sigma2_0), lam0]])
        res = minimize(negloglik, x0, method='Nelder-Mead', options={'maxiter': 40000, 'xatol': 1e-10, 'fatol': 1e-10})
        if best is None or res.fun < best.fun:
            best = res
    return -best.fun, best.x[:k]


def dea_correlation_check(sub, extra_col):
    sub = sub.dropna(subset=[extra_col])
    sub = sub[sub[extra_col] > 0].copy()
    r2, r3 = [], []
    for yr, g in sub.groupby('year'):
        if len(g) < 5:
            continue
        X2 = g[['cows_heads', 'fodder_sown_area_ha']].values.astype(float)
        X3 = g[['cows_heads', 'fodder_sown_area_ha', extra_col]].values.astype(float)
        Y = g[['milk_production_tons']].values.astype(float)
        r2.append(pd.Series(dea_input_oriented(X2, Y, vrs=True), index=g.index))
        r3.append(pd.Series(dea_input_oriented(X3, Y, vrs=True), index=g.index))
    sub['vrs2'] = pd.concat(r2); sub['vrs3'] = pd.concat(r3)
    r, p = pearsonr(sub['vrs2'], sub['vrs3'])
    return sub, r, p


def run():
    aug = data_prep.build_augmented_panel()

    # ---- Table 7: Labor ----
    print('=== Table 7: Labor-augmented DEA ===')
    lab_sub, r_lab, p_lab = dea_correlation_check(aug, 'total_labor_hours')
    print(f'N={len(lab_sub)}, mean2={lab_sub["vrs2"].mean():.3f}, mean3={lab_sub["vrs3"].mean():.3f}, r={r_lab:.3f} (p={p_lab:.2e})')
    pd.DataFrame({'measure': ['N', 'mean_2input', 'mean_3input', 'correlation', 'p_value'],
                  'value': [len(lab_sub), lab_sub['vrs2'].mean(), lab_sub['vrs3'].mean(), r_lab, p_lab]
                  }).to_csv(os.path.join(OUT_DIR, 'table8_labor_robustness.csv'), index=False)

    # Waldman (1982) diagnostic: 3-input SFA has positively-skewed OLS residual
    y = np.log(lab_sub['milk_production_tons'].values)
    X3 = np.column_stack([np.ones(len(lab_sub)), np.log(lab_sub['cows_heads'].values),
                          np.log(lab_sub['fodder_sown_area_ha'].values), np.log(lab_sub['total_labor_hours'].values)])
    beta_ols = np.linalg.lstsq(X3, y, rcond=None)[0]
    resid_skew = skew(y - X3 @ beta_ols)
    print(f'Waldman diagnostic: 3-input OLS residual skewness = {resid_skew:.3f} (positive = problematic for SFA)')

    # ---- SI S12: Capital ----
    print('\n=== SI S12: Capital-augmented ===')
    cap_sub, r_cap, p_cap = dea_correlation_check(aug, 'capital_stock_thsd_manat')
    print(f'N={len(cap_sub)}, mean2={cap_sub["vrs2"].mean():.3f}, mean3={cap_sub["vrs3"].mean():.3f}, r={r_cap:.3f} (p={p_cap:.2e})')

    y = np.log(cap_sub['milk_production_tons'].values)
    lx1 = np.log(cap_sub['cows_heads'].values); lx2 = np.log(cap_sub['fodder_sown_area_ha'].values)
    lx4 = np.log(cap_sub['capital_stock_thsd_manat'].values)
    n = len(y)
    logL2, _ = sfa_loglik(y, np.column_stack([np.ones(n), lx1, lx2]))
    logL3, beta3 = sfa_loglik(y, np.column_stack([np.ones(n), lx1, lx2, lx4]))
    LR_cap = 2 * (logL3 - logL2)
    p_lr_cap = 1 - chi2.cdf(LR_cap, 1)
    print(f'SFA LR test (capital): chi2(1)={LR_cap:.3f}, p={p_lr_cap:.4f}, capital coef={beta3[3]:.4f}')
    pd.DataFrame({'measure': ['N', 'mean_2input', 'mean_3input', 'correlation', 'p_value_corr', 'LR_stat', 'LR_p', 'capital_coef'],
                  'value': [len(cap_sub), cap_sub['vrs2'].mean(), cap_sub['vrs3'].mean(), r_cap, p_cap, LR_cap, p_lr_cap, beta3[3]]
                  }).to_csv(os.path.join(OUT_DIR, 's12_capital_robustness.csv'), index=False)

    # ---- SI S13: Land ----
    print('\n=== SI S13: Land-augmented ===')
    land_sub, r_land, p_land = dea_correlation_check(aug, 'land_use_ha')
    print(f'N={len(land_sub)}, mean2={land_sub["vrs2"].mean():.3f}, mean3={land_sub["vrs3"].mean():.3f}, r={r_land:.3f} (p={p_land:.2e})')

    y = np.log(land_sub['milk_production_tons'].values)
    lx1 = np.log(land_sub['cows_heads'].values); lx2 = np.log(land_sub['fodder_sown_area_ha'].values)
    lx5 = np.log(land_sub['land_use_ha'].values)
    n = len(y)
    logL2b, _ = sfa_loglik(y, np.column_stack([np.ones(n), lx1, lx2]))
    logL3b, beta3b = sfa_loglik(y, np.column_stack([np.ones(n), lx1, lx2, lx5]))
    LR_land = 2 * (logL3b - logL2b)
    p_lr_land = 1 - chi2.cdf(LR_land, 1)
    print(f'SFA LR test (land): chi2(1)={LR_land:.3f}, p={p_lr_land:.4f}, land coef={beta3b[3]:.4f}')
    pd.DataFrame({'measure': ['N', 'mean_2input', 'mean_3input', 'correlation', 'p_value_corr', 'LR_stat', 'LR_p', 'land_coef'],
                  'value': [len(land_sub), land_sub['vrs2'].mean(), land_sub['vrs3'].mean(), r_land, p_land, LR_land, p_lr_land, beta3b[3]]
                  }).to_csv(os.path.join(OUT_DIR, 's13_land_robustness.csv'), index=False)

    # ---- Joint capital + land ----
    print('\n=== Joint capital + land ===')
    full = aug.dropna(subset=['capital_stock_thsd_manat', 'land_use_ha'])
    full = full[(full['capital_stock_thsd_manat'] > 0) & (full['land_use_ha'] > 0)]
    y = np.log(full['milk_production_tons'].values)
    lx1 = np.log(full['cows_heads'].values); lx2 = np.log(full['fodder_sown_area_ha'].values)
    lx4 = np.log(full['capital_stock_thsd_manat'].values); lx5 = np.log(full['land_use_ha'].values)
    n = len(y)
    logL2c, _ = sfa_loglik(y, np.column_stack([np.ones(n), lx1, lx2]))
    logL4c, beta4c = sfa_loglik(y, np.column_stack([np.ones(n), lx1, lx2, lx4, lx5]))
    LR_joint = 2 * (logL4c - logL2c)
    p_joint = 1 - chi2.cdf(LR_joint, 2)
    print(f'N={n}, LR joint: chi2(2)={LR_joint:.3f}, p={p_joint:.4f}, betas={beta4c}')
    pd.DataFrame({'measure': ['N', 'districts', 'LR_stat', 'LR_p', 'intercept', 'cattle_coef', 'fodder_coef', 'capital_coef', 'land_coef'],
                  'value': [n, full['region'].nunique(), LR_joint, p_joint, *beta4c]
                  }).to_csv(os.path.join(OUT_DIR, 's13_joint_capital_land.csv'), index=False)


if __name__ == '__main__':
    run()
