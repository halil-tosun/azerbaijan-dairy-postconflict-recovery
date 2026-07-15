"""
00_data_prep.py
================
Builds every analytical panel used in the paper, starting from the raw
official statistics files in ../data/.

DATA SOURCES
------------
- ../data/original_panel/azerbaijan_dairy_panel_WIDE_districts.csv
      District-level panel (milk, cattle, fodder, cost price, profitability,
      labour-hours/centner). Pre-compiled from State Statistical Committee
      of the Republic of Azerbaijan district-level yearbook tables.
- ../data/original_panel/azerbaijan_dairy_panel_WIDE_economic_regions.csv
      Economic-region-level panel (13 regions x 2000-2024, complete), used
      for the Malmquist decomposition (Table 5).
- ../data/official_stats_raw/002_2en.xls
      "Land in ownership and(or) use, ha" -- enterprise subsector, by
      district/region, State Statistical Committee table 2.2.
- ../data/official_stats_raw/003_5en.xls, 003_6en.xls
      "Main production funds of agriculture purposes" (capital stock) --
      private-farm subsector, by district/region, tables 3.5 and 3.6.
- ../data/official_stats_raw/002_53en.xls
      "Labour expenditure for per centner of milk, person-hour" --
      enterprise subsector, by district/region, table 2.53.

DOCUMENTED DATA-CLEANING DECISIONS (see paper Section 3.1 and SI S12)
----------------------------------------------------------------------
1. Sample = 67 administrative districts (cities excluded).
2. Aghdara district (1 usable observation, 2024) is RETAINED in the main
   production panel (N=1,523) but EXCLUDED from district-ranking exercises
   (SI Table S8), which require >=20 years of data per district.
3. Milk cost price is coalesced from enterprise + private-farm series;
   five extreme values (Tukey far-out fence) are treated as missing in the
   main analysis: Gobustan 2013, Khachmaz 2016, Khizi 2015, Jalilabad 2015,
   Shabran 2019. All three treatments (retain all / exclude one / exclude
   five) are reproduced (SI Table S9).
4. Capital stock (private-farm subsector) and land-use (enterprise
   subsector) are DISTINCT SUBSECTOR MEASURES, not all-category totals;
   see SI Section S12 for the full data-scope discussion.
5. Post-conflict / event-study treated districts: Aghdam, Fuzuli, Gubadli,
   Jabrayil, Kalbajar, Khojavand, Lachin, Shusha, Zangilan (9 districts;
   Aghdara excluded, single observation only).
"""
import pandas as pd
import numpy as np
import csv
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'output')
os.makedirs(OUT_DIR, exist_ok=True)

AGHDARA = 'Aghdara district'
MIN_YEARS_FOR_RANKING = 20

COST_OUTLIERS = [
    ('Gobustan district', 2013), ('Khachmaz district', 2016),
    ('Khizi district', 2015), ('Jalilabad district', 2015),
    ('Shabran district', 2019),
]

POST_CONFLICT_DISTRICTS_ALL = [
    'Aghdam district', 'Aghdara district', 'Fuzuli district', 'Gubadli district',
    'Jabrayil district', 'Kalbajar district', 'Khojavand district', 'Lachin district',
    'Shusha district', 'Zangilan district',
]
# Event-study treated group excludes Aghdara (single observation only)
EVENT_STUDY_TREATED = [d for d in POST_CONFLICT_DISTRICTS_ALL if d != AGHDARA]


def load_district_panel(exclude_aghdara=False, outlier_treatment='exclude_five'):
    """Cleaned district-level panel (one row per district-year)."""
    raw = pd.read_csv(os.path.join(DATA_DIR, 'original_panel', 'azerbaijan_dairy_panel_WIDE_districts.csv'))
    d = raw[raw['region_type'] == 'district'].copy()
    if exclude_aghdara:
        d = d[d['region'] != AGHDARA]

    d['cost_milk'] = d['cost_price_per_centner_milk_enterprises_manat'].combine_first(
        d['cost_price_per_centner_milk_private_manat'])
    d['profitability'] = d['profitability_milk_enterprises_pct'].combine_first(
        d['profitability_milk_private_pct'])

    if outlier_treatment == 'exclude_five':
        mask = pd.Series(False, index=d.index)
        for region, year in COST_OUTLIERS:
            mask |= (d['region'] == region) & (d['year'] == year)
        d.loc[mask, 'cost_milk'] = np.nan
    elif outlier_treatment == 'exclude_gobustan_only':
        mask = (d['region'] == 'Gobustan district') & (d['year'] == 2013)
        d.loc[mask, 'cost_milk'] = np.nan
    elif outlier_treatment == 'retain_all':
        pass
    else:
        raise ValueError(outlier_treatment)
    return d


def dea_analysis_sample(d):
    """District-years with complete output/input data (N=1,523 production panel)."""
    return d.dropna(subset=['milk_production_tons', 'cows_heads', 'fodder_sown_area_ha']).copy()


def load_region_panel():
    """Economic-region-level panel (13 regions x 25 years, complete)."""
    df = pd.read_csv(os.path.join(DATA_DIR, 'original_panel', 'azerbaijan_dairy_panel_WIDE_economic_regions.csv'))
    return df.dropna(subset=['milk_production_tons', 'cows_heads', 'fodder_sown_area_ha']).copy()


def _parse_regional_xls_csv(path, value_name):
    """Parse a State Statistical Committee regional table (already exported to CSV via LibreOffice)."""
    with open(path, encoding='utf-8') as f:
        rows = list(csv.reader(f))
    header_idx = None
    for i, row in enumerate(rows):
        if any(c.strip().strip('"') == '2005' for c in row):
            header_idx = i
            break
    years = [int(y) for y in rows[header_idx][2:] if y.strip().replace('.', '').isdigit()]
    records = []
    for row in rows[header_idx + 1:]:
        if len(row) < 3:
            continue
        region = row[1].strip()
        if not region:
            continue
        for yr, v in zip(years, row[2:2 + len(years)]):
            v = v.strip()
            if v in ('-', '...', ''):
                continue
            v = v.replace(',', '.')
            try:
                val = float(v)
            except ValueError:
                continue
            records.append({'region_raw': region, 'year': yr, value_name: val})
    return pd.DataFrame(records)


def load_labor_data():
    """
    District-level labour expenditure per centner of milk (person-hours),
    converted to total labour-hours. Source: table 2.53 (enterprise subsector).
    NOTE: requires the .xls to be pre-converted to CSV (see README /
    convert_xls_to_csv.sh) since this environment's pandas has no xlrd.
    """
    path_csv = os.path.join(DATA_DIR, 'official_stats_raw', '002_53en.csv')
    labor = _parse_regional_xls_csv(path_csv, 'labour_hours_per_centner')
    labor['region_raw'] = labor['region_raw'].str.strip()
    return labor.rename(columns={'region_raw': 'region'})


def load_capital_data():
    """District-level agricultural capital stock (thousand AZN), private-farm subsector. Table 3.5."""
    path_csv = os.path.join(DATA_DIR, 'official_stats_raw', '003_5en.csv')
    cap = _parse_regional_xls_csv(path_csv, 'capital_stock_thsd_manat')
    cap['region_raw'] = cap['region_raw'].str.strip()
    return cap.rename(columns={'region_raw': 'region'})


def load_land_data():
    """District-level land in ownership/use (ha), enterprise subsector. Table 2.2."""
    path_csv = os.path.join(DATA_DIR, 'official_stats_raw', '002_2en.csv')
    land = _parse_regional_xls_csv(path_csv, 'land_use_ha')
    land['region_raw'] = land['region_raw'].str.strip()
    return land.rename(columns={'region_raw': 'region'})


def build_augmented_panel():
    """Main DEA/SFA sample merged with labor, capital, and land (for Section 3.7 robustness checks)."""
    d = dea_analysis_sample(load_district_panel(exclude_aghdara=False))
    labor = load_labor_data()
    capital = load_capital_data()
    land = load_land_data()

    d = d.merge(labor, on=['region', 'year'], how='left')
    d['total_labor_hours'] = d['labour_hours_per_centner'] * d['milk_production_tons'] * 10
    d = d.merge(capital, on=['region', 'year'], how='left')
    d = d.merge(land, on=['region', 'year'], how='left')
    return d


if __name__ == '__main__':
    d = dea_analysis_sample(load_district_panel(exclude_aghdara=False))
    print('Districts:', d['region'].nunique(), '| N (production panel):', len(d))
    d.to_csv(os.path.join(OUT_DIR, 'analysis_panel_districts.csv'), index=False)

    regions = load_region_panel()
    print('Economic regions:', regions['region'].nunique(), '| obs:', len(regions))
    regions.to_csv(os.path.join(OUT_DIR, 'analysis_panel_regions.csv'), index=False)

    aug = build_augmented_panel()
    print('Augmented panel (labor/capital/land merged):', len(aug))
    print('  N with labor:', aug['total_labor_hours'].notna().sum())
    print('  N with capital:', aug['capital_stock_thsd_manat'].notna().sum())
    print('  N with land:', aug['land_use_ha'].notna().sum())
    aug.to_csv(os.path.join(OUT_DIR, 'analysis_panel_augmented.csv'), index=False)

    # Table 1 (paper): Description of variables used in the empirical analysis
    table1 = pd.DataFrame([
        ('Milk production', 'Annual milk output', 'Tons', 'Output', 'SSC'),
        ('Dairy cattle stock', 'Dairy cattle inventory', 'Head', 'Input', 'SSC'),
        ('Fodder area', 'Fodder crop sown area', 'Hectares', 'Input', 'SSC'),
        ('Milk production cost', 'Average production cost', 'AZN per ton', 'Efficiency determinant', 'SSC'),
        ('Milk profitability', 'Profitability of milk production', 'Percent', 'Efficiency determinant', 'SSC'),
        ('Labor input*', 'Total labor hours', 'Hours', 'Robustness analysis', 'SSC'),
        ('Agricultural capital*', 'Main production funds', 'AZN', 'Robustness analysis', 'SSC'),
        ('Agricultural land*', 'Land in ownership/use', 'Hectares', 'Robustness analysis', 'SSC'),
    ], columns=['Variable', 'Definition', 'Unit', 'Role', 'Source'])
    table1.to_csv(os.path.join(OUT_DIR, 'table1_variable_description.csv'), index=False)
    print('\nTable 1 (variable description) written.')
