"""Step 3 - hot pixels.

A hot pixel is a detector pixel that reads high in most spectra of a session, always at the same wavelength.
For each species: every spectrum minus its 7-pixel running median gives a residual; a pixel spikes when the residual
is above 5 robust standard deviations (1.4826 x MAD). A pixel that spikes in more than 30 % of the spectra of the
species is hot. Hot pixels and their immediate neighbours are replaced by linear interpolation.
"""
import numpy as np
from scipy.signal import medfilt

from .config import HOT_PIXEL_FRACTION, HOT_PIXEL_MEDIAN, HOT_PIXEL_Z


def detect_hot_pixels(wl, spectra):
    """Return the boolean mask of pixels to replace and the wavelengths (nm) of the hot pixels."""
    resid = spectra - np.array([medfilt(x, HOT_PIXEL_MEDIAN) for x in spectra])
    sigma = np.median(np.abs(resid), 1, keepdims=True) * 1.4826 + 1e-9
    frac = (resid / sigma > HOT_PIXEL_Z).mean(0)
    hot = np.where(frac > HOT_PIXEL_FRACTION)[0]
    mask = np.zeros(wl.size, bool)
    for i in hot:
        mask[max(0, i - 1):i + 2] = True
    return mask, [round(float(wl[i]), 2) for i in hot]


def interpolate_pixels(wl, x, mask):
    x = x.copy()
    x[mask] = np.interp(wl[mask], wl[~mask], x[~mask])
    return x
