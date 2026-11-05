"""Compare the Python descriptors (02_colour_maps.py) with pavo (03_pavo_descriptors.R).

Usage:
    python 04_check_pavo.py <data folder>

Prints the largest absolute difference of B1, B2, B3, H1, S6 and S8 over every geometry, and the number of
geometries where H1 differs. Writes <data folder>/pavo_check.csv.
"""
import os
import sys

import pandas as pd

data = sys.argv[1] if len(sys.argv) == 2 else sys.exit(__doc__)
py = pd.concat([pd.read_csv(os.path.join(data, 'metrics_by_geometry.csv')).assign(v='main'),
                pd.read_csv(os.path.join(data, 'metrics_previous_sweeps.csv')).assign(v='previous_sweep')])
pv = pd.read_csv(os.path.join(data, 'pavo_descriptors.csv')).rename(columns={'variant': 'v'})
m = py.merge(pv, on=['species_id', 'v', 'I_deg', 'C_deg'], suffixes=('_python', '_pavo'))
assert len(m) == len(py) == len(pv), (len(m), len(py), len(pv))
rows = []
for k in ['B1', 'B2', 'B3', 'H1', 'S6', 'S8']:
    d = (m[f'{k}_python'] - m[f'{k}_pavo']).abs()
    rows.append(dict(descriptor=k, geometries=len(m), max_abs_difference=d.max(), n_above_0_001=int((d > 1e-3).sum())))
r = pd.DataFrame(rows)
r.to_csv(os.path.join(data, 'pavo_check.csv'), index=False)
print(r.to_string(index=False))
