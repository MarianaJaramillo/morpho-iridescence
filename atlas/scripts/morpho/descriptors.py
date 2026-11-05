"""Montgomerie (2006) colour descriptors, with the definitions of pavo::summary.rspec (Maia et al. 2013, 2019).

Computed on the cleaned spectrum over 400-700 nm, with negative values set to zero
(the same as pavo::procspec(fixneg = 'zero')). R is in % of the white reference.

  B1  total brightness   sum of R                      (%.nm, on the 1 nm grid)
  B2  mean brightness    mean of R                     (%)
  B3  intensity          max R                         (%)
  H1  peak wavelength    wavelength of max R           (nm)
  S6  contrast           max R - min R                 (percentage points)
  S8  chroma             (max R - min R) / B2          (dimensionless)

The same values are produced by scripts/03_pavo_descriptors.R; scripts/04_check_pavo.py compares both.
"""
import numpy as np

from .config import WL_VIS


def montgomerie(grid, spectra):
    v = (grid >= WL_VIS[0]) & (grid <= WL_VIS[1])
    wl = grid[v]
    R = np.clip(spectra[:, v], 0, None)
    B1, B2, B3, Rmin = R.sum(1), R.mean(1), R.max(1), R.min(1)
    with np.errstate(divide='ignore', invalid='ignore'):
        S8 = np.where(B2 > 0, (B3 - Rmin) / B2, np.nan)
    return dict(B1=B1, B2=B2, B3=B3, H1=wl[R.argmax(1)].astype(float), S6=B3 - Rmin, S8=S8)
