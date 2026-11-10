"""Colour maps and descriptors of every geometry, from the cleaned spectra.

Usage:
    python 02_colour_maps.py <data folder>

Reads cleaned_spectra/*.csv and spectra_index.csv (written by 01_clean_spectra.py) and computes, per geometry:
    * Montgomerie / pavo descriptors B1, B2, B3, H1, S6, S8          morpho/descriptors.py
    * colour under CIE D65, 2-degree observer:
        perceived colour (common exposure) and CIE hue (chromaticity only)   morpho/cie_colour.py
    * measurement checks (row drift, negative baseline, suspect labels)    morpho/qa.py

Outputs in <data folder>:
    metrics_by_geometry.csv           main datasets (32 species), one row per (species, I, C)
    metrics_previous_sweeps.csv       the earlier sweeps of the 5 remeasured species
    dataset_summary.csv               one row per dataset: coverage, noise, row drift, specular peak, tilt
    coverage_grid.csv                 every I x C position: measured, blind_zone or missing
    measurement_checks_log.csv        row drift, suspect labels and missing geometries found
    parameters.json                   every parameter used, including the colour exposure
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from morpho import changelog, config
from morpho.cie_colour import cie_hue, linear_srgb, perceived_colour, xyz
from morpho.config import C_GRID, I_GRID, PREVIOUS, in_blind_zone
from morpho.descriptors import montgomerie
from morpho.qa import measurement_flags, tilt_estimate
from morpho.resample import GRID


def main(data):
    index = pd.read_csv(os.path.join(data, 'spectra_index.csv'))
    hot = pd.read_csv(os.path.join(data, 'hot_pixels.csv')).set_index(['species_id', 'variant']).hot_pixels_nm
    metas, rgbs, Ys, summary, coverage = [], [], [], [], []
    for (sid, variant), meta in index.groupby(['species_id', 'variant'], sort=False):
        suffix = '' if variant != PREVIOUS else '__previous_sweep'
        W = pd.read_csv(os.path.join(data, 'cleaned_spectra', f'{sid}{suffix}.csv'))
        meta = meta.reset_index(drop=True)
        assert (W.I_deg.values == meta.I_deg.values).all() and (W.C_deg.values == meta.C_deg.values).all()
        CL = W.drop(columns=['I_deg', 'C_deg']).values
        for k, v in montgomerie(GRID, CL).items():
            meta[k] = v
        meta = measurement_flags(meta, GRID, CL, sid, variant)
        XYZ = xyz(CL)
        rgbs.append(linear_srgb(XYZ)); Ys.append(XYZ[:, 1])
        tilt, tilt_iqr = tilt_estimate(meta)
        have = set(zip(meta.I_deg, meta.C_deg))
        missing = []
        for I in I_GRID:
            for C in C_GRID:
                status = 'measured' if (I, C) in have else ('blind_zone' if in_blind_zone(I, C) else 'missing')
                coverage.append(dict(species_id=sid, species=meta.species.iloc[0], variant=variant, I_deg=int(I), C_deg=int(C), status=status))
                if status == 'missing':
                    missing.append(int(I))
        if missing:
            byI = pd.Series(missing).value_counts().sort_index().to_dict()
            changelog.log(sid, variant, '-', 'missing', f'{len(missing)} geometries outside the blind zone were not measured. By I: {byI}')
        band = meta[((meta.C_deg + meta.I_deg).abs() <= 10) & (meta.I_deg.abs() <= 50)]
        summary.append(dict(species_id=sid, species=meta.species.iloc[0], variant=variant, n_spectra=len(meta), n_missing=len(missing),
                            hot_pixels_nm=hot.get((sid, variant), ''), median_noise_pp=meta.noise_sigma_pp.median(),
                            fraction_row_drift=meta.flag_row_drift.mean(), max_row_offset_pp=meta.row_offset_pp.abs().max(),
                            fraction_negative_baseline=meta.flag_negative_baseline.mean(), n_suspect_labels=int(meta.flag_suspect_label.sum()),
                            specular_B3_median_pct=band.B3.median(), specular_H1_median_nm=band.H1.median(),
                            tilt_deg=tilt, tilt_iqr_deg=tilt_iqr))
        metas.append(meta)
    M = pd.concat(metas, ignore_index=True)
    rgb = np.vstack(rgbs)
    exposure = float(np.percentile(np.concatenate(Ys), config.EXPOSURE_PERCENTILE))
    M['perceived_colour_hex'] = perceived_colour(rgb, exposure)
    M['cie_hue_hex'] = cie_hue(rgb)
    cols = ['species_id', 'species', 'variant', 'I_deg', 'C_deg', 'B1', 'B2', 'B3', 'H1', 'S6', 'S8',
            'perceived_colour_hex', 'cie_hue_hex', 'source_file', 'source_column', 'original_header',
            'integration_time_ms', 'n_averages', 'n_spikes_removed', 'noise_sigma_pp', 'row_offset_pp',
            'flag_row_drift', 'flag_negative_baseline', 'flag_above_white_std', 'flag_suspect_label']
    M = M[cols].sort_values(['species', 'variant', 'I_deg', 'C_deg'])
    prev = M.variant == PREVIOUS
    M[~prev].round(4).to_csv(os.path.join(data, 'metrics_by_geometry.csv'), index=False)
    M[prev].round(4).to_csv(os.path.join(data, 'metrics_previous_sweeps.csv'), index=False)
    pd.DataFrame(summary).round(3).to_csv(os.path.join(data, 'dataset_summary.csv'), index=False)
    pd.DataFrame(coverage).to_csv(os.path.join(data, 'coverage_grid.csv'), index=False)
    changelog.to_frame().to_csv(os.path.join(data, 'measurement_checks_log.csv'), index=False)
    params = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in vars(config).items()
              if k.isupper() and k not in ('SPECIES', 'REMEASURED')}
    params['REMEASURED'] = sorted(config.REMEASURED)
    params['exposure_Y'] = exposure
    json.dump(params, open(os.path.join(data, 'parameters.json'), 'w'), indent=2)
    print(f'geometries: {(~prev).sum()} main + {prev.sum()} previous sweeps | exposure Y = {exposure:.4f}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
