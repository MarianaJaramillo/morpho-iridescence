"""Build the Morpho Angular Atlas (HTML viewer) from the data folder.

Usage:
    python 07_build_viewer.py <data folder> <viewer output folder>

Writes:
    <viewer>/web/morpho_atlas.html + web/spectra/*.txt   for a web server (spectra load on demand)
    <viewer>/Morpho_Angular_Atlas_offline.html           one self-contained file: open it by double-clicking

Template: scripts/viewer/atlas_template.html. The data are inserted where the template says /*__DATA__*/null.
Spectra are stored as base64 int16 (R x 10, every second wavelength): cleaned 380-780 nm at 2 nm, then raw.
"""
import base64
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

if len(sys.argv) != 3:
    sys.exit(__doc__)
DATA, OUT = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
WEB = os.path.join(OUT, 'web')
shutil.rmtree(WEB, ignore_errors=True)
os.makedirs(os.path.join(WEB, 'spectra'))

M = pd.concat([pd.read_csv(os.path.join(DATA, 'metrics_by_geometry.csv')),
               pd.read_csv(os.path.join(DATA, 'metrics_previous_sweeps.csv'))], ignore_index=True)
order = []
for sid in sorted(M.species_id.unique()):
    for v in ['single', 'remeasurement', 'previous_sweep']:
        if ((M.species_id == sid) & (M.variant == v)).any():
            order.append((sid, v))

datasets, files = [], []
for sid, v in order:
    g = M[(M.species_id == sid) & (M.variant == v)].sort_values(['I_deg', 'C_deg'])
    name = g.species.iloc[0]
    label = {'single': name, 'remeasurement': name + ' (remeasured)', 'previous_sweep': name + ' (previous sweep)'}[v]
    srcfiles = sorted(g.source_file.unique())
    suffix = '__previous_sweep' if v == 'previous_sweep' else ''
    CL = pd.read_csv(os.path.join(DATA, 'cleaned_spectra', f'{sid}{suffix}.csv'))
    RW = pd.read_csv(os.path.join(DATA, 'raw_spectra', f'{sid}{suffix}.csv'))
    assert (CL.I_deg.values == g.I_deg.values).all() and (RW.C_deg.values == g.C_deg.values).all()
    raw_wl = RW.columns[2:].astype(float).values
    rows = []
    nz = lambda x, d=1: None if pd.isna(x) else round(float(x), d)
    for _, r in g.iterrows():
        flags = int(r.flag_row_drift) << 1 | int(r.flag_negative_baseline) << 2 | int(r.flag_above_white_std) << 3 | int(r.flag_suspect_label) << 4
        # row layout read by the template: 0 I, 1 C, 2 perceived colour, 3 CIE hue, 4 B2, 5 H1, 6 S8, 7 B3, 8 flags,
        # 9 row offset, 10-12 unused, 13 source file index, 14 original header, 15 unused, 16 S6
        rows.append([int(r.I_deg), int(r.C_deg), r.perceived_colour_hex, r.cie_hue_hex, round(float(r.B2), 2), int(r.H1),
                     nz(r.S8, 2), round(float(r.B3), 1), flags, round(float(r.row_offset_pp), 1), None, None, None,
                     srcfiles.index(r.source_file), str(r.original_header).replace('.RFL8', ''), None, round(float(r.S6), 1)])
    c = CL.values[:, 2:][:, ::2]; rr = RW.values[:, 2:][:, ::2]
    blob = np.clip(np.round(np.concatenate([c, rr], 1) * 10), -32000, 32000).astype('<i2').tobytes()
    fn = f'spectra/{sid}__{v}.txt'
    open(os.path.join(WEB, fn), 'w').write(base64.b64encode(blob).decode())
    files.append(fn)
    datasets.append(dict(id=f'{sid}__{v}', sid=sid, variant=v, name=name, label=label, file=fn, n=len(g), srcfiles=srcfiles, rows=rows))

cov = pd.read_csv(os.path.join(DATA, 'coverage_grid.csv'))
missing = {f'{sid}__{v}': [[int(a), int(b)] for a, b in zip(g.I_deg, g.C_deg)]
           for (sid, v), g in cov[cov.status == 'missing'].groupby(['species_id', 'variant'])}
gruson = json.load(open(os.path.join(DATA, 'gruson', 'gruson_fits.json')))
align = json.load(open(os.path.join(DATA, 'gruson', 'gruson_alignment.json')))
payload = dict(gruson=gruson, gruson_align=align, datasets=datasets, missing=missing,
               cleanWl=[380, 2, c.shape[1]], rawWl=[round(float(x), 2) for x in raw_wl[::2]], nClean=c.shape[1], nRaw=rr.shape[1])
tpl = open(os.path.join(HERE, 'viewer', 'atlas_template.html')).read()
html = tpl.replace('/*__DATA__*/null', json.dumps(payload, separators=(',', ':')))
open(os.path.join(WEB, 'morpho_atlas.html'), 'w').write(html)

# offline version: the spectra are embedded as <script type="text/plain" id="spec-<dataset>"> before the main script
blocks = ''.join(f'<script type="text/plain" id="spec-{os.path.basename(f)[:-4]}">{open(os.path.join(WEB, f)).read()}</script>\n' for f in files)
page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n</head>\n<body>\n'
        + blocks + html.replace('<meta charset="utf-8">\n', '', 1) + '\n</body>\n</html>\n')
open(os.path.join(OUT, 'Morpho_Angular_Atlas_offline.html'), 'w').write(page)
print(f'datasets {len(datasets)} | web page {len(html) / 1e6:.2f} MB + {len(files)} spectra files | offline page {len(page) / 1e6:.1f} MB')
