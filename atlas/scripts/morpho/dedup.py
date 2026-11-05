"""Step 2 - deduplication.

  * Exact copies: two columns with the same spectrum (same hash of the values rounded to 4 decimals), inside a file
    or across the positive and negative files. One is kept.
  * Same geometry measured twice with different content: the column from the file whose sign matches I is kept
    (negative file for I <= 0, positive file for I > 0).

A final check verifies that every kept spectrum is still the one of its original column.
"""
import hashlib

import numpy as np

from .changelog import log


def spectrum_hash(x):
    return hashlib.md5(np.nan_to_num(x).round(4).tobytes()).hexdigest()


def deduplicate(meta, spectra, species_id, variant):
    meta = meta.copy()
    meta['keep'] = True
    # exact copies
    for _, g in meta[meta.duplicated('hash', keep=False)].groupby('hash'):
        meta.loc[g.index[1:], 'keep'] = False
    n_copies = int((~meta.keep).sum())
    if n_copies:
        per = meta[~meta.keep].groupby('source_file').size().to_dict()
        log(species_id, variant, ';'.join(per), 'identical_copies',
            f'{n_copies} columns are exact copies of another spectrum of the same species and were dropped: {per}')
    # same geometry, different content
    k = meta[meta.keep]
    for (I, C), g in k[k.duplicated(['I', 'C'], keep=False)].groupby(['I', 'C']):
        pref = g[g.from_negative_file == (I <= 0)]
        keep = pref.index[0] if len(pref) else g.index[0]
        meta.loc[g.index.difference([keep]), 'keep'] = False
        r = np.corrcoef(spectra[g.index][:, 400:900])[0, 1]
        log(species_id, variant, ';'.join(g.source_file.unique()), 'duplicate_geometry',
            f'I={I:g} C={C:g} measured {len(g)} times (r={r:.2f}); kept "{meta.original_header[keep]}" from {meta.source_file[keep]}.')
    idx = meta.index[meta.keep].values
    meta = meta.loc[idx].drop(columns='keep').reset_index(drop=True)
    spectra = spectra[idx]
    assert all(spectrum_hash(spectra[i]) == meta.hash[i] for i in range(len(meta))), 'deduplication misaligned the spectra'
    return meta, spectra
