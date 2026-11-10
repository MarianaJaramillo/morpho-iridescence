"""Cleaning of the raw angular reflectance spectra (steps 1-6).

Usage:
    python 01_clean_spectra.py [--only <species_id>[,<species_id>...]] <Spectra folder> <remeasurements folder> <corrections folder> <output data folder>

    --only    process just those species and merge them into the tables already in <output data folder>
              (cleaned_spectra/, raw_spectra/, spectra_index.csv, hot_pixels.csv, cleaning_log.csv)

    <Spectra folder>          the 64 xlsx files of Spectra.zip plus the 2 of M. hecuba (one positive and one negative file per species)
    <remeasurements folder>   the 10 xlsx files remeasured for deidamia, eugenia, lympharis, menelaus and rhetenor
    <corrections folder>      the corrected negative files of cisseis and amphitryon (replace the ones in Spectra)

For each species, and for both sweeps of the 5 remeasured species:
    1. header repair             morpho/headers.py
    2. deduplication             morpho/dedup.py
    3. hot pixels                morpho/hot_pixels.py
    4. Hampel filter             morpho/hampel.py
    5. Savitzky-Golay smoothing  morpho/smoothing.py
    6. resampling 380-780 nm     morpho/resample.py

Outputs in <output data folder>:
    cleaned_spectra/<species_id>.csv           I_deg, C_deg, 380 ... 780 (R in %, 1 nm)
    raw_spectra/<species_id>.csv               I_deg, C_deg, raw wavelengths in 380-780 nm (R in %, angles repaired, nothing else changed)
    (previous sweeps end in __previous_sweep.csv)
    spectra_index.csv                          one row per spectrum: geometry, original file, column and header, acquisition settings
    hot_pixels.csv                             hot pixels found in each dataset
    cleaning_log.csv                           every correction applied, with file and column of origin
The raw xlsx files are read only; they are never modified.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from morpho import changelog
from morpho.config import PREVIOUS, REMEASURED, REMEASUREMENT, SINGLE, SPECIES, WL_OUT, WL_WORK
from morpho.dedup import deduplicate, spectrum_hash
from morpho.hampel import hampel, noise_sigma
from morpho.headers import assign_off_grid, repair_file_headers
from morpho.hot_pixels import detect_hot_pixels, interpolate_pixels
from morpho.io_raw import file_sign, files_for, read_avantes_xlsx
from morpho.resample import GRID, resample
from morpho.smoothing import savitzky_golay


def load_dataset(species_id, pattern, folder, variant, corrections):
    """Steps 1 and 2: read all files of one dataset, repair headers, remove duplicates."""
    recs, spec, WL = [], [], None
    for f in files_for(pattern, folder, species_id, variant, corrections):
        fname = os.path.basename(f)
        wl, X, names, it, av = read_avantes_xlsx(f)
        WL = wl if WL is None else WL
        assert np.allclose(wl, WL)
        for j, ((I, C), n) in enumerate(zip(repair_file_headers(names, species_id, variant, fname), names)):
            recs.append(dict(species_id=species_id, variant=variant, I=I, C=C, source_file=fname, source_column=j + 1,
                             original_header=str(n), integration_time_ms=it[j], n_averages=av[j],
                             hash=spectrum_hash(X[j]), from_negative_file=file_sign(f) == 'negative'))
            spec.append(X[j])
    meta = pd.DataFrame(recs)
    spec = np.array(spec)
    meta = assign_off_grid(meta, species_id, variant)
    meta, spec = deduplicate(meta, spec, species_id, variant)
    return meta, spec, WL


def clean_dataset(meta, spec, wl, species_id, variant):
    """Steps 3-6 on the working range 350-820 nm, then resampling to 380-780 nm."""
    w_mask = (wl >= WL_WORK[0]) & (wl <= WL_WORK[1])
    w = wl[w_mask]
    Xw = spec[:, w_mask]
    bad, hot_nm = detect_hot_pixels(w, Xw)
    changelog.log(species_id, variant, 'all', 'hot_pixels',
                  f'Pixels with a systematic spike: {hot_nm} nm (plus their immediate neighbours); interpolated.')
    cleaned, extra = [], []
    for x in Xw:
        x = interpolate_pixels(w, x, bad)              # step 3
        x, n_spikes = hampel(x)                        # step 4
        sm = savitzky_golay(x)                         # step 5
        cleaned.append(resample(w, sm))                # step 6
        extra.append(dict(n_spikes_removed=n_spikes, noise_sigma_pp=noise_sigma(w, x)))
    meta = pd.concat([meta.reset_index(drop=True), pd.DataFrame(extra)], axis=1)
    return meta, np.array(cleaned), hot_nm


def merge_table(path, new, only):
    """Replace the rows of the species in `only` in an existing table, keeping the order of config.SPECIES."""
    old = pd.read_csv(path)
    pos = {sid: k for k, (sid, _, _) in enumerate(SPECIES)}
    merged = pd.concat([old[~old.species_id.isin(only)], new], ignore_index=True)
    return merged.iloc[merged.species_id.map(pos).argsort(kind='stable')]


def main(spectra_dir, remeasure_dir, corrections_dir, out, only=None):
    os.makedirs(os.path.join(out, 'cleaned_spectra'), exist_ok=True)
    os.makedirs(os.path.join(out, 'raw_spectra'), exist_ok=True)
    jobs = []
    for sid, name, pat in SPECIES:
        if only and sid not in only:
            continue
        if sid in REMEASURED:
            jobs += [(sid, name, pat, remeasure_dir, REMEASUREMENT, None), (sid, name, pat, spectra_dir, PREVIOUS, corrections_dir)]
        else:
            jobs.append((sid, name, pat, spectra_dir, SINGLE, corrections_dir))
    index, hot = [], []
    for sid, name, pat, folder, variant, corr in jobs:
        meta, spec, wl = load_dataset(sid, pat, folder, variant, corr)
        meta, CL, hot_nm = clean_dataset(meta, spec, wl, sid, variant)
        meta['species'] = name
        meta = meta.rename(columns={'I': 'I_deg', 'C': 'C_deg'})
        meta[['I_deg', 'C_deg']] = meta[['I_deg', 'C_deg']].astype(int)
        order = meta.sort_values(['I_deg', 'C_deg']).index.values
        meta, CL, spec = meta.loc[order].reset_index(drop=True), CL[order], spec[order]
        suffix = '' if variant != PREVIOUS else '__previous_sweep'
        W = pd.DataFrame(CL.round(4), columns=[str(int(x)) for x in GRID])
        W.insert(0, 'C_deg', meta.C_deg.values); W.insert(0, 'I_deg', meta.I_deg.values)
        W.to_csv(os.path.join(out, 'cleaned_spectra', f'{sid}{suffix}.csv'), index=False)
        rv = (wl >= WL_OUT[0]) & (wl <= WL_OUT[1])
        Rw = pd.DataFrame(spec[:, rv].round(4), columns=[f'{x:.2f}' for x in wl[rv]])
        Rw.insert(0, 'C_deg', meta.C_deg.values); Rw.insert(0, 'I_deg', meta.I_deg.values)
        Rw.to_csv(os.path.join(out, 'raw_spectra', f'{sid}{suffix}.csv'), index=False)
        index.append(meta); hot.append(dict(species_id=sid, variant=variant, hot_pixels_nm=';'.join(f'{h:.2f}' for h in hot_nm)))
        print(f'{sid:13s} {variant:15s} spectra={len(meta):4d} hot pixels={hot_nm}', flush=True)
    cols = ['species_id', 'species', 'variant', 'I_deg', 'C_deg', 'source_file', 'source_column', 'original_header',
            'integration_time_ms', 'n_averages', 'n_spikes_removed', 'noise_sigma_pp']
    idx = pd.concat(index, ignore_index=True)[cols]
    idx['noise_sigma_pp'] = idx.noise_sigma_pp.round(4)
    hot, log = pd.DataFrame(hot), changelog.to_frame()
    if only:
        idx = merge_table(os.path.join(out, 'spectra_index.csv'), idx, only)
        hot = merge_table(os.path.join(out, 'hot_pixels.csv'), hot, only)
        log = merge_table(os.path.join(out, 'cleaning_log.csv'), log, only)
    idx.to_csv(os.path.join(out, 'spectra_index.csv'), index=False)
    hot.to_csv(os.path.join(out, 'hot_pixels.csv'), index=False)
    log.to_csv(os.path.join(out, 'cleaning_log.csv'), index=False)
    print('spectra:', len(idx), '| datasets:', len(jobs))


if __name__ == '__main__':
    args = sys.argv[1:]
    only = None
    if args[:1] == ['--only']:
        if len(args) < 2:
            sys.exit(__doc__)
        only, args = set(args[1].split(',')), args[2:]
    if len(args) != 4:
        sys.exit(__doc__)
    main(*args, only=only)
