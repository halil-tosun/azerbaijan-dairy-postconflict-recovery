"""
run_all.py
==========
Runs the full replication pipeline in order and writes all output
tables to ../output/.

Expected runtime: ~30-40 minutes total, dominated by the two bootstrap
procedures:
  - 01_dea_bootstrap.py:      ~15 min (B=200, 25 annual cross-sections)
  - 04_malmquist_bootstrap.py: ~15-20 min (B=300, 24 year-pairs x 13 regions)
Both print progress after each year / year-pair. 04 saves incrementally
and can be re-run to resume if interrupted (pass start_idx/end_idx).

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

    print('\n=== 01: DEA + Simar-Wilson bootstrap (Table 3, Figure 2) ===')
    m01 = _load('01_dea_bootstrap')
    m01.run()

    print('\n=== 02: SFA models -- translog PREFERRED (Table 4, SI S1-S2) ===')
    m02 = _load('02_sfa_models')
    m02.run()

    print('\n=== 03: Second-stage regression + robustness (Table 6, SI S3/S5/S6/S9) ===')
    m03 = _load('03_second_stage')
    m03.run()

    print('\n=== 04: Bootstrap Malmquist (Table 5, Figures 3-5) ===')
    m04 = _load('04_malmquist_bootstrap')
    m04.run()

    print('\n=== 05: Event-study (Table 7, SI S10) ===')
    m05 = _load('05_event_study')
    m05.run()

    print('\n=== 06: Labor/capital/land robustness (Table 8, SI S12-S13) ===')
    m06 = _load('06_labor_capital_land_robustness')
    m06.run()

    print('\n=== 07: Diagnostics (SI S4/S7/S8) ===')
    m07 = _load('07_robustness_diagnostics')
    m07.run()

    print('\n=== 08a: Figure 1 (analytical framework diagram, 600 DPI) ===')
    m08a = _load('08a_make_figure1_framework')
    m08a.run()

    print('\n=== 08: Figures 2-6 (600 DPI, from output/ CSVs) ===')
    m08 = _load('08_make_figures')
    m08.run()

    print(f'\nAll done in {time.time()-t0:.0f} seconds. See ../output/ for all result files and ../figures/ for all figures.')
