"""
02_sfa_models.py
=================
Stochastic Frontier Analysis. Translog is the PAPER'S PREFERRED MODEL
(Table 3), selected over Cobb-Douglas via a likelihood-ratio test
(chi2(3)=15.67, p=0.0013). Cobb-Douglas (SI Table S1) and a panel
model with a linear time trend (SI Table S2) are reported as
robustness checks.

Model: Aigner, Lovell & Schmidt (1977) normal-half-normal composed
error frontier, estimated by maximum likelihood.
"""
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, chi2
import importlib.util
import os

spec = importlib.util.spec_from_file_location("data_prep", os.path.join(os.path.dirname(__file__), "00_data_prep.py"))
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)
OUT_DIR = data_prep.OUT_DIR


def _mle_halfnormal(y, X, lam_starts=(0.5, 1.0, 1.5, 2.0, 3.0, 5.0)):
    n, k = X.shape

    def negloglik(params):
        beta = params[:k]
        sigma2 = np.exp(params[k])
        lam = params[k + 1]
        sigma = np.sqrt(sigma2)
        eps = y - X @ beta
        z = eps / sigma
        return -np.sum(np.log(2) - np.log(sigma) + norm.logpdf(z) + norm.logcdf(-z * lam))

    beta0 = np.linalg.lstsq(X, y, rcond=None)[0]
    sigma2_0 = np.var(y - X @ beta0)
    best = None
    for lam0 in lam_starts:
        x0 = np.concatenate([beta0, [np.log(sigma2_0), lam0]])
        res = minimize(negloglik, x0, method='Nelder-Mead',
                        options={'maxiter': 40000, 'xatol': 1e-10, 'fatol': 1e-10})
        if best is None or res.fun < best.fun:
            best = res
    beta = best.x[:k]
    sigma2 = np.exp(best.x[k])
    lam = best.x[k + 1]
    sigma_v2 = sigma2 / (1 + lam ** 2)
    sigma_u2 = sigma2 - sigma_v2
    gamma = sigma_u2 / sigma2
    logL = -best.fun

    eps = y - X @ beta
    sigma = np.sqrt(sigma2)
    sigma_star = np.sqrt(sigma_u2 * sigma_v2 / sigma2)
    mu_star = -eps * sigma_u2 / sigma2
    Eu = mu_star + sigma_star * norm.pdf(mu_star / sigma_star) / norm.cdf(mu_star / sigma_star)
    TE = np.exp(-Eu)
    return dict(beta=beta, sigma2=sigma2, lam=lam, sigma_u2=sigma_u2, sigma_v2=sigma_v2,
                gamma=gamma, logL=logL, TE=TE)


def translog_sfa(d):
    """Table 3 (PREFERRED / MAIN MODEL)."""
    y = np.log(d['milk_production_tons'].values)
    lx1 = np.log(d['cows_heads'].values)
    lx2 = np.log(d['fodder_sown_area_ha'].values)
    X = np.column_stack([np.ones(len(y)), lx1, lx2, 0.5 * lx1 ** 2, 0.5 * lx2 ** 2, lx1 * lx2])
    return _mle_halfnormal(y, X)


def cobb_douglas_sfa(d):
    """SI Table S1 (functional-form robustness)."""
    y = np.log(d['milk_production_tons'].values)
    X = np.column_stack([np.ones(len(y)), np.log(d['cows_heads'].values), np.log(d['fodder_sown_area_ha'].values)])
    return _mle_halfnormal(y, X)


def panel_trend_sfa(d):
    """SI Table S2 (autonomous technical-change robustness)."""
    y = np.log(d['milk_production_tons'].values)
    lx1 = np.log(d['cows_heads'].values)
    lx2 = np.log(d['fodder_sown_area_ha'].values)
    trend = d['year'].values - d['year'].values.min()
    X = np.column_stack([np.ones(len(y)), lx1, lx2, trend])
    return _mle_halfnormal(y, X)


def run():
    d = data_prep.dea_analysis_sample(data_prep.load_district_panel(exclude_aghdara=False))

    print('--- Table 3: Translog SFA (PREFERRED MODEL) ---')
    tl = translog_sfa(d)
    print('beta:', tl['beta'], 'logL:', tl['logL'])
    print('TE mean/sd/min/max:', tl['TE'].mean(), tl['TE'].std(), tl['TE'].min(), tl['TE'].max())

    print('\n--- SI Table S1: Cobb-Douglas SFA (robustness) ---')
    cd = cobb_douglas_sfa(d)
    print('beta:', cd['beta'], 'logL:', cd['logL'])

    LR = 2 * (tl['logL'] - cd['logL'])
    pval = 1 - chi2.cdf(LR, df=3)
    print(f'\nLR test (translog vs Cobb-Douglas): chi2(3) = {LR:.3f}, p = {pval:.4f}')

    print('\n--- SI Table S2: Panel SFA with linear time trend (robustness) ---')
    pt = panel_trend_sfa(d)
    print('beta:', pt['beta'], 'logL:', pt['logL'])

    d = d.copy()
    d['sfa_te_translog'] = tl['TE']
    d.to_csv(os.path.join(OUT_DIR, 'table4_sfa_district_te.csv'), index=False)

    pd.DataFrame({
        'model': ['Translog (Table 3, PREFERRED)', 'Cobb-Douglas (S1)', 'Panel+trend (S2)'],
        'logL': [tl['logL'], cd['logL'], pt['logL']],
        'sigma2': [tl['sigma2'], cd['sigma2'], pt['sigma2']],
        'lambda': [tl['lam'], cd['lam'], pt['lam']],
        'gamma': [tl['gamma'], cd['gamma'], pt['gamma']],
    }).to_csv(os.path.join(OUT_DIR, 'table4_s1_s2_sfa_summary.csv'), index=False)

    return tl, cd, pt


if __name__ == '__main__':
    run()
