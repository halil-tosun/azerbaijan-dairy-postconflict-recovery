"""
08_make_figures.py
====================
Generates all main-text and Supporting Information figures from the CSV
outputs produced by scripts 01, 04, and 09. Run after run_all.py (or after
the relevant individual scripts).

This script's plotting parameters (figure size, colors, 300 DPI) are the
exact ones used to produce the figures embedded in the submitted manuscript
and Supporting Information; it is the single authoritative figure-generation
script for this replication package.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import t as tdist
import os

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'output')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans'})


def figure2_annual_efficiency():
    """Figure 2 (main text): annual mean bootstrap-corrected DEA efficiency."""
    fig2 = pd.read_csv(os.path.join(OUT_DIR, 'figure2_annual_efficiency.csv'))
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(fig2['year'], fig2['vrs_bc'], marker='o', markersize=3, linewidth=1.5, color='#1b4965')
    ax.set_xlabel('Year')
    ax.set_ylabel('Mean bootstrap-corrected DEA technical efficiency')
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure2_Annual_Efficiency.png'), dpi=300)
    plt.close()
    print('Figure 2 saved as Figure2_Annual_Efficiency.png')


def figure3_regional_tfp_timeseries():
    """Figure 3 (main text): regional cumulative Malmquist TFP indices."""
    fig3 = pd.read_csv(os.path.join(OUT_DIR, 'figure3_regional_tfp_timeseries.csv'))
    fig, ax = plt.subplots(figsize=(7, 5))
    regions = fig3['region'].unique()
    cmap = plt.cm.tab20(np.linspace(0, 1, len(regions)))
    for i, r in enumerate(regions):
        sub = fig3[fig3['region'] == r].sort_values('year')
        label = r.replace(' economic region', '').replace(' Autonomous Republic', ' A.R.')
        ax.plot(sub['year'], sub['cum_TFP'], linewidth=1.3, color=cmap[i], label=label)
    ax.axhline(1.0, color='gray', linestyle=':', linewidth=1)
    ax.set_xlabel('Year')
    ax.set_ylabel('Cumulative Malmquist TFP index (2000 = 1)')
    ax.legend(fontsize=6.5, ncol=2, loc='upper left', framealpha=0.9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure3_Regional_TFP_TimeSeries.png'), dpi=300)
    plt.close()
    print('Figure 3 saved as Figure3_Regional_TFP_TimeSeries.png')


def figure4_event_study():
    """Figure 4 (main text): event-study coefficients, control-frontier
    (primary) and pooled-frontier (secondary) dependent variables."""
    cf = pd.read_csv(os.path.join(OUT_DIR, 'figure6_event_study_control_frontier.csv')).sort_values('event_time')
    pf = pd.read_csv(os.path.join(OUT_DIR, 's10_event_study_pooled_frontier.csv')).sort_values('event_time')

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), sharey=True)
    for ax, df, title in [(axes[0], cf, 'Control-frontier score (primary)'),
                           (axes[1], pf, 'Pooled-frontier score (secondary)')]:
        crit = tdist.ppf(0.975, 1418)
        ci = crit * df['se_clustered']
        ax.errorbar(df['event_time'], df['coef'], yerr=ci, fmt='o', markersize=4,
                     capsize=3, color='#1b4965', ecolor='#5fa8d3', linewidth=1)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.axvline(-0.5, color='red', linestyle='--', linewidth=1, alpha=0.6)
        ax.set_xlabel('Event time (years since 2020 reintegration)')
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel('Coefficient rel. to k=\u22121\n(95% CI, clustered SE)')
    fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.13, wspace=0.08)
    plt.savefig(os.path.join(FIG_DIR, 'Figure4_EventStudy.png'), dpi=300)
    plt.close()
    print('Figure 4 saved as Figure4_EventStudy.png')


def figures1_2_regional_decomposition():
    """Figure S1 (Supporting Information): regional TC/TEC decomposition."""
    t5 = pd.read_csv(os.path.join(OUT_DIR, 'table5_summary.csv')).sort_values('cum_TFP_pct', ascending=True)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    y = range(len(t5))
    labels = [r.replace(' economic region', '').replace(' Autonomous Republic', ' A.R.') for r in t5['region']]
    ax.barh(y, t5['cum_TC_pct'], color='#1b4965', label='Technological change (TC)', height=0.4, align='edge')
    ax.barh([p + 0.4 for p in y], t5['cum_TEC_pct'], color='#bc4749',
            label='Technical efficiency change (TEC)', height=0.4, align='edge')
    ax.set_yticks([p + 0.4 for p in y])
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Cumulative change, 2000\u20132024 (%)')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'FigureS1_Regional_Decomposition.png'), dpi=300)
    plt.close()
    print('Figure S1 (SI) saved as FigureS1_Regional_Decomposition.png')


def figures2_technology_efficiency_matrix():
    """Figure S2 (Supporting Information): regional technology-efficiency matrix."""
    t5 = pd.read_csv(os.path.join(OUT_DIR, 'table5_summary.csv'))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(t5['cum_TC_pct'], t5['cum_TEC_pct'], s=40, color='#1b4965', zorder=3)
    for _, row in t5.iterrows():
        label = row['region'].replace(' economic region', '').replace(' Autonomous Republic', ' A.R.')
        ax.annotate(label, (row['cum_TC_pct'], row['cum_TEC_pct']), fontsize=6.5,
                    xytext=(4, 3), textcoords='offset points')
    ax.axhline(0, color='gray', linewidth=0.8)
    ax.axvline(0, color='gray', linewidth=0.8)
    ax.set_xlabel('Cumulative technological change, TC (%)')
    ax.set_ylabel('Cumulative technical efficiency change, TEC (%)')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'FigureS2_Technology_Efficiency_Matrix.png'), dpi=300)
    plt.close()
    print('Figure S2 (SI) saved as FigureS2_Technology_Efficiency_Matrix.png')


def run():
    figure2_annual_efficiency()
    figure3_regional_tfp_timeseries()
    figure4_event_study()
    figures1_2_regional_decomposition()
    figures2_technology_efficiency_matrix()


if __name__ == '__main__':
    run()
