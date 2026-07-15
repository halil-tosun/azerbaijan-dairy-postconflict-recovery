"""
08_make_figures.py
====================
Generates all 5 figures from the CSV outputs produced by scripts 01-05.
Run after run_all.py (or after the relevant individual scripts).
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'output')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)


def figure2_dea():
    fig1 = pd.read_csv(os.path.join(OUT_DIR, 'figure2_annual_efficiency.csv'))
    plt.figure(figsize=(9, 5.5))
    plt.plot(fig1['year'], fig1['vrs_bc'], marker='o', color='#1f77b4', linewidth=2, markersize=6)
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('Mean bootstrap bias-corrected DEA efficiency', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure2_DEA_Annual_Efficiency.png'), dpi=600)
    plt.close()
    print('Figure 1 saved')


def figure3_malmquist_ts():
    ts = pd.read_csv(os.path.join(OUT_DIR, 'figure3_regional_tfp_timeseries.csv'))
    pivot = ts.pivot(index='year', columns='region', values='cum_TFP')
    pivot.columns = [c.replace(' economic region', '').replace('Nakhchivan Autonomous Republic', 'Nakhchivan A.R.') for c in pivot.columns]
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)
    plt.figure(figsize=(10, 6.2))
    colors = plt.cm.tab20(range(len(pivot.columns)))
    for i, col in enumerate(pivot.columns):
        plt.plot(pivot.index, pivot[col], label=col, linewidth=1.6, color=colors[i])
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('Cumulative TFP Index (2000=1)', fontsize=12)
    plt.legend(loc='upper left', fontsize=7.5, ncol=2, frameon=True)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure3_Regional_TFP_TimeSeries.png'), dpi=600)
    plt.close()
    print('Figure 2 saved')


def figure4_decomp():
    t5 = pd.read_csv(os.path.join(OUT_DIR, 'table5_summary.csv')).sort_values('cum_TFP_pct', ascending=False)
    labels = [r.replace(' economic region', '').replace('Nakhchivan Autonomous Republic', 'Nakhchivan A.R.') for r in t5['region']]
    x = np.arange(len(labels)); width = 0.27
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(x - width, t5['cum_TFP_pct'], width, label='TFP (cumulative, %)', color='#1a3a5c')
    ax.bar(x, t5['cum_TC_pct'], width, label='Technological Change, TC (%)', color='#2e86ab')
    ax.bar(x + width, t5['cum_TEC_pct'], width, label='Technical Efficiency Change, TEC (%)', color='#d9622b')
    for i, v in enumerate(t5['cum_TFP_pct']):
        ax.text(i - width, v + (8 if v >= 0 else -14), f"{v:+.0f}", ha='center', fontsize=8, fontweight='bold', color='#1a3a5c')
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=40, ha='right', fontsize=9)
    ax.set_ylabel('Cumulative change, 2000-2024 (%)', fontsize=12)
    ax.axhline(0, color='black', linewidth=0.8)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.25, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure4_Regional_Decomposition.png'), dpi=600)
    plt.close()
    print('Figure 3 saved')


def figure5_matrix():
    t5 = pd.read_csv(os.path.join(OUT_DIR, 'table5_summary.csv'))
    labels = [r.replace(' economic region', '').replace('Nakhchivan Autonomous Republic', 'Nakhchivan A.R.') for r in t5['region']]
    x = t5['cum_TEC_pct'].values; y = t5['cum_TC_pct'].values
    size = np.abs(t5['cum_TFP_pct'].values); color = t5['cum_TFP_pct'].values
    fig, ax = plt.subplots(figsize=(9.5, 8))
    sc = ax.scatter(x, y, s=size * 3.2 + 80, c=color, cmap='RdYlBu_r', edgecolors='black', linewidths=0.7, alpha=0.9, vmin=-100, vmax=300)
    for i, lab in enumerate(labels):
        ax.annotate(lab, (x[i], y[i]), fontsize=8, fontweight='bold', ha='center', va='center')
    lims = [min(x.min(), y.min()) - 20, max(x.max(), y.max()) + 20]
    ax.plot(lims, lims, '--', color='gray', linewidth=1.2)
    ax.axhline(0, color='black', linewidth=0.8); ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Technical Efficiency Change, TEC (cumulative %, 2000-2024)', fontsize=11)
    ax.set_ylabel('Technological Change, TC (cumulative %, 2000-2024)', fontsize=11)
    cb = plt.colorbar(sc, ax=ax); cb.set_label('Cumulative TFP growth (%)', fontsize=10)
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure5_Technology_Efficiency_Matrix.png'), dpi=600)
    plt.close()
    print('Figure 4 saved')


def figure6_eventstudy():
    ev = pd.read_csv(os.path.join(OUT_DIR, 'figure6_event_study.csv'))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = ev['event_time'].values; y = ev['coef'].values; ci = 1.96 * ev['se'].values
    ax.errorbar(x, y, yerr=ci, fmt='o-', color='#1a3a5c', ecolor='#7fa8c9', capsize=3, linewidth=1.8, markersize=5)
    ax.axhline(0, color='gray', linewidth=0.8)
    ax.axvline(-0.5, color='red', linestyle='--', linewidth=1.2, label='2020 reintegration')
    ax.set_xlabel('Years relative to reintegration (2020 = 0)', fontsize=12)
    ax.set_ylabel('Efficiency gap vs. control districts\n(relative to event time = -1)', fontsize=12)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure6_EventStudy.png'), dpi=600)
    plt.close()
    print('Figure 5 saved')


def run():
    figure2_dea(); figure3_malmquist_ts(); figure4_decomp(); figure5_matrix(); figure6_eventstudy()


if __name__ == '__main__':
    run()
