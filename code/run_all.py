"""
run_all.py
==========
Runs the full replication pipeline in order and writes all output
tables to ../output/.

Expected runtime: ~60-75 minutes total, dominated by the bootstrap
procedures:
  - 01_dea_bootstrap.py:                ~30-33 min (B=200, 25 annual cross-sections)
  - 04_malmquist_bootstrap.py:          ~15-20 min (B=300, 24 year-pairs x 13 regions)
  - 09_control_frontier_event_study.py: ~7-9 min (control-only-frontier bootstrap DEA)
  - 10_gpt_review_robustness.py:        ~9-11 min (dominated by its own Aghdara-excluded
                                          re-run of the control-frontier construction)
All bootstrap steps print progress after each year / year-pair. 04 saves
incrementally and can be re-run to resume if interrupted (pass start_idx/end_idx).

Script 10 formalizes four additional robustness/diagnostic analyses added to
the paper in response to an external methodological review (SI Tables S15-S18):
DEA-SFA rank correlation, a second-stage specification with year fixed effects,
control-frontier bootstrap/truncation diagnostics, and an Aghdara-exclusion +
wild-cluster-bootstrap check on the primary causal estimate. It must run after
01, 02, 03, and 09.

Run individual numbered scripts directly to regenerate only one part.
"""
import importlib.util
import os
import time

HERE = os.path.dirname(__file__)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


if __name__ == '__main__':
    t0 = time.time()
    print('=== 00: data preparation ===')
    _load('00_data_prep')

    print('\n=== 01: DEA + Simar-Wilson bootstrap (Table 2, Figure 2) ===')
    m01 = _load('01_dea_bootstrap')
    m01.run()

    print('\n=== 02: SFA models -- translog PREFERRED (Table 3, SI S1-S2) ===')
    m02 = _load('02_sfa_models')
    m02.run()

    print('\n=== 03: Second-stage regression + robustness (Table 5, SI S3.1/S5/S6/S9) ===')
    m03 = _load('03_second_stage')
    m03.run()

    print('\n=== 04: Bootstrap Malmquist (Table 4, Figure 3, SI Figures S1-S2) ===')
    m04 = _load('04_malmquist_bootstrap')
    m04.run()

    print('\n=== 05: Event-study (pooled-frontier, legacy/diagnostic only -- superseded by 09) ===')
    m05 = _load('05_event_study')
    m05.run()

    print('\n=== 09: Control-frontier event-study, clustered SE, fractional logit, MDE (Tables 6-8, Figure 4) ===')
    m09 = _load('09_control_frontier_event_study')
    m09.run()

    print('\n=== 06: Labor/capital/land robustness (Table 9, SI S10-S12) ===')
    m06 = _load('06_labor_capital_land_robustness')
    m06.run()

    print('\n=== 07: Diagnostics (SI S4/S7/S8) ===')
    m07 = _load('07_robustness_diagnostics')
    m07.run()

    print('\n=== 10: Additional robustness from external review (SI S15-S18) ===')
    m10 = _load('10_gpt_review_robustness')
    m10.run_s15_rank_correlation()
    m10.run_s16_year_fe_second_stage()
    m10.run_s17_control_frontier_diagnostics()
    m10.run_s18b_wild_cluster_bootstrap()
    m10.run_s18a_aghdara_sensitivity()  # slowest (~9-11 min); run last

    print('\n=== 08: Figures 2-4 + SI Figures S1-S2 (600 DPI, from output/ CSVs) ===')
    m08 = _load('08_make_figures')
    m08.run()

    print(f'\nAll done in {time.time()-t0:.0f} seconds. See ../output/ for all result files and ../figures/ for all figures.')
    print('\nNote: 08a_make_figure1_framework.py is NOT part of the active pipeline.')
    print('It produces an analytical-framework diagram that is not referenced')
    print('anywhere in the current manuscript body (Methods, Results, or SI) and')
    print('is retained in the repository only as an optional, unused legacy script.')
