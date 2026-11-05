"""Place the Gruson measurement points on the goniometer maps.

Usage:
    python 06_gruson_alignment.py <data folder>

The Gruson sessions use a different angle convention. A Gruson point (I_g, C_g) is drawn on the map at
    I = b + s/2,  C = b - s/2,   with  b = -(I_g - C_g)/2 + d,  s = -(I_g + C_g)
i.e. the map orientation is (-I_g, C_g); this orientation correlates best with the map brightness (median r 0.80
against 0.37 for the direct one). d is a shift along the bisector, chosen per dataset (steps of 0.25 degrees,
-20 to 20) to maximise the correlation between the Gruson brightness points and B2 interpolated on the map.
helenor results are linked only to Morpho helenor (UNAL).

Writes <data folder>/gruson/gruson_alignment.json and .csv: per dataset, d (degrees), r at d = 0 and r at d.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

data = sys.argv[1] if len(sys.argv) == 2 else sys.exit(__doc__)
G = json.load(open(os.path.join(data, 'gruson', 'gruson_fits.json')))
M = pd.concat([pd.read_csv(os.path.join(data, 'metrics_by_geometry.csv')),
               pd.read_csv(os.path.join(data, 'metrics_previous_sweeps.csv'))])
Is, Cs = np.arange(-80, 81, 10), np.arange(-90, 91, 5)


def b2_interpolator(g):
    A = np.full((17, 37), np.nan)
    for I, C, b in zip(g.I_deg, g.C_deg, g.B2):
        A[(int(I) + 80) // 10, (int(C) + 90) // 5] = b
    for i in range(17):                     # fill the blind zone along each row
        ok = ~np.isnan(A[i])
        if ok.sum() > 1:
            A[i] = np.interp(Cs, Cs[ok], A[i][ok])
    return RegularGridInterpolator((Is, Cs), A, bounds_error=False, fill_value=np.nan)


out = {}
for (sid, variant), g in M.groupby(['species_id', 'variant']):
    key = sid if sid in G else ('helenor' if sid == 'helenor_unal' else None)
    if not key:
        continue
    f = b2_interpolator(g)
    P = [(p[2], p[0] + p[1], p[3]) for sweep in G[key]['B'] for p in sweep['pts']]
    x = np.array([p[0] for p in P]); s = np.array([p[1] for p in P]); y = np.array([p[2] for p in P])

    def r_at(d):
        b = -x + d; sp = -s
        v = f(np.c_[b + sp / 2, b - sp / 2]); ok = ~np.isnan(v)
        return np.corrcoef(v[ok], y[ok])[0, 1] if ok.sum() > 4 else np.nan

    ds = np.arange(-20, 20.25, 0.25)
    rs = np.array([r_at(d) for d in ds]); i = np.nanargmax(rs)
    out[f'{sid}__{variant}'] = dict(d=float(ds[i]), r0=round(float(r_at(0)), 3), r=round(float(rs[i]), 3))
json.dump(out, open(os.path.join(data, 'gruson', 'gruson_alignment.json'), 'w'), indent=1)
T = pd.DataFrame(out).T.rename_axis('dataset').rename(columns={'d': 'shift_d_deg', 'r0': 'r_at_d0', 'r': 'r_at_shift'})
T.to_csv(os.path.join(data, 'gruson', 'gruson_alignment.csv'))
print(T.to_string()); print('median r at d=0:', T.r_at_d0.median(), '| median |d|:', T.shift_d_deg.abs().median())
