"""Reading the Avantes exports (.RFL8 saved as .xlsx) and finding the files of each species.

Layout of every raw file (read without modifying it):
    row 0   spectrum names, e.g. 'M. absoloni -10I 20C.RFL8' (one column per geometry)
    row 2   integration time (ms)
    row 3   number of averages
    row 6+  column 0 = wavelength (nm), columns 1.. = reflectance (% of the white reference)
"""
import glob
import os
import re

import numpy as np
import pandas as pd

from .changelog import log


def read_avantes_xlsx(path):
    df = pd.read_excel(path, header=None)
    wavelength = df.iloc[6:, 0].astype(float).values
    spectra = df.iloc[6:, 1:].apply(pd.to_numeric, errors='coerce').values.T
    names = df.iloc[0, 1:].tolist()
    integration_ms = df.iloc[2, 1:].tolist()
    n_averages = df.iloc[3, 1:].tolist()
    return wavelength, spectra, names, integration_ms, n_averages


def file_sign(path):
    """'negative' for files of negative incidence angles (the raw names also use the misspelling 'negaitivo')."""
    b = os.path.basename(path).lower()
    return 'negative' if ('negativ' in b or 'negaitiv' in b) else 'positive'


def _matching(pattern, folder):
    files = [f for f in sorted(glob.glob(os.path.join(folder, '*.xlsx')))
             if re.search(pattern, os.path.basename(f).replace('_', ' '))]
    return [f for f in files if not (pattern == 'M. helenor p' and 'unal' in f)]


def files_for(pattern, folder, species_id, variant, corrections_folder=None):
    """Raw files of one species in one folder.

    A file in `corrections_folder` replaces the file of the same species and the same sign
    (used for the negative files of cisseis and amphitryon, which in Spectra.zip were copies of the positive file).
    """
    files = _matching(pattern, folder)
    if corrections_folder:
        for c in _matching(pattern, corrections_folder):
            old = [f for f in files if file_sign(f) == file_sign(c)]
            files = [f for f in files if f not in old] + [c]
            log(species_id, variant, os.path.basename(c), 'file_replaced',
                f'{os.path.basename(c)} is used instead of {", ".join(os.path.basename(o) for o in old)} '
                f'(the original file was a copy of the positive file).')
    return sorted(files)
