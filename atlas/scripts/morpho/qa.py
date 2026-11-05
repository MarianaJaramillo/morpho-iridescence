"""Measurement checks. These are quality-control indicators of the measurement, not colour descriptors.

  row_offset_pp          median B2 of the incidence row away from specular (|C + I| >= 40) minus the mean of the
                         neighbouring rows (percentage points)
  flag_row_drift         |row_offset_pp| > 3 pp: the whole row is shifted up or down (reference drift)
  flag_negative_baseline median of the spectrum over 400-700 nm below -2 % (dark-reference offset)
  flag_above_white_std   B3 > 100 %: brighter than the diffuse white standard (expected near specular)
  flag_suspect_label     cell away from specular (|C + I| > 15) with B2 >= 20 % and more than 5 times the median of its
                         neighbours: probably a spectrum labelled with the wrong geometry. Not moved.
  tilt_deg               (C_peak + I) / 2, median over rows |I| <= 60, where C_peak is the parabolic peak of B2 in the row
"""
import numpy as np

from .changelog import log
from .config import NEGATIVE_BASELINE_PCT, ROW_DRIFT_MAX_PP, SUSPECT_FACTOR, WL_VIS


def measurement_flags(meta, grid, spectra, species_id, variant):
    v = (grid >= WL_VIS[0]) & (grid <= WL_VIS[1])
    meta['flag_negative_baseline'] = np.median(spectra[:, v], 1) < NEGATIVE_BASELINE_PCT
    meta['flag_above_white_std'] = meta.B3 > 100
    off = meta[(meta.C_deg + meta.I_deg).abs() >= 40].groupby('I_deg').B2.median()
    nb = ((off.shift(1) + off.shift(-1)) / 2).fillna(off.shift(1)).fillna(off.shift(-1))
    meta['row_offset_pp'] = meta.I_deg.map(off - nb)
    meta['flag_row_drift'] = meta.row_offset_pp.abs() > ROW_DRIFT_MAX_PP
    for I, g in meta[meta.flag_row_drift].groupby('I_deg'):
        log(species_id, variant, g.source_file.iloc[0], 'row_drift',
            f'I={I:g}: reflectance away from specular is shifted {g.row_offset_pp.iloc[0]:+.1f} pp from the neighbouring rows. Not corrected.')
    t = meta.set_index(['I_deg', 'C_deg']).B2
    sus = []
    for (I, C), b in t.items():
        if abs(C + I) <= 15 or b < 20:
            sus.append(False); continue
        nbv = [t.get(k) for k in [(I, C - 5), (I, C + 5), (I - 10, C), (I + 10, C), (I, C - 10), (I, C + 10)]]
        nbv = [x for x in nbv if x is not None and np.isfinite(x)]
        sus.append(len(nbv) >= 2 and b > SUSPECT_FACTOR * max(np.median(nbv), 2.0))
    meta['flag_suspect_label'] = sus
    for _, r in meta[meta.flag_suspect_label].iterrows():
        near = [t.get(k) for k in [(r.I_deg, r.C_deg - 5), (r.I_deg, r.C_deg + 5), (r.I_deg - 10, r.C_deg), (r.I_deg + 10, r.C_deg)] if t.get(k) is not None]
        log(species_id, variant, r.source_file, 'suspect_label',
            f'"{r.original_header}": B2 {r.B2:.0f} % against ~{np.median(near):.0f} % in its neighbours; probably belongs to another geometry. Not moved.',
            r.source_column)
    return meta


def tilt_estimate(meta):
    pts = []
    for I, g in meta[meta.I_deg.abs() <= 60].groupby('I_deg'):
        g = g.sort_values('C_deg'); y = g.B2.values; k = y.argmax(); c = g.C_deg.values
        if 0 < k < len(g) - 1 and c[k + 1] - c[k - 1] == 10:
            den = y[k - 1] - 2 * y[k] + y[k + 1]
            off = 0.5 * (y[k - 1] - y[k + 1]) / den if den else 0
            pts.append((I, c[k] + off * 5))
    if len(pts) < 4:
        return np.nan, np.nan
    d = np.array([(c + i) / 2 for i, c in pts])
    return float(np.median(d)), float(np.subtract(*np.percentile(d, [75, 25])))
