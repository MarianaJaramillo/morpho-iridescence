"""Spectral cleaning: smoothing, artefact masking, white-reference normalisation, descriptors."""
import numpy as np
from scipy.signal import savgol_filter


def smooth(y, window=151, polyorder=3):
    y = y.copy()
    bad = np.isnan(y)
    if bad.any():
        y[bad] = np.interp(np.flatnonzero(bad), np.flatnonzero(~bad), y[~bad])
    return savgol_filter(y, window, polyorder)


def mask_band(wl, y, band):
    out = y.copy()
    m = (wl >= band[0]) & (wl <= band[1])
    if m.any() and (~m).any():
        out[m] = np.interp(wl[m], wl[~m], out[~m])
    return out


def nearest_reference(I, C, ref_geometries):
    """Reference closest in (I, C) by L1 distance; ties go to the lowest (I, C)."""
    return min(ref_geometries, key=lambda g: (abs(I - g[0]) + abs(C - g[1]), g))


def hue_peak_centroid(wl, y, frac=0.10):
    """Intensity-weighted centroid of the region within `frac` of the peak height (above the minimum)."""
    ipk = np.argmax(y)
    thr = y[ipk] - frac * (y[ipk] - y.min())
    sel = y >= thr
    if sel.sum() < 2:
        return wl[ipk]
    return np.sum(wl[sel] * y[sel]) / np.sum(y[sel])


def hue_edge(wl, y, plateau=(420, 470), floor=(620, 670)):
    """Wavelength of the half-fall between the blue plateau and the red floor (for spectra without a peak)."""
    top = np.median(y[(wl >= plateau[0]) & (wl <= plateau[1])])
    bot = np.median(y[(wl >= floor[0]) & (wl <= floor[1])])
    mid = (top + bot) / 2
    i0 = np.argmax(wl >= plateau[1])
    idx = np.where(y[i0:] < mid)[0]
    if len(idx) == 0 or idx[0] == 0:
        return np.nan
    j = i0 + idx[0]
    return np.interp(mid, [y[j], y[j - 1]], [wl[j], wl[j - 1]])
