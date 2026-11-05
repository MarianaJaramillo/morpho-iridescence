"""Step 4 - Hampel filter (random spikes).

Within each spectrum, a point is a spike when it differs from the 7-pixel running median by more than
6 robust standard deviations (1.4826 x MAD of the residuals of that spectrum). Spikes are replaced by the median.
"""
import numpy as np
from scipy.signal import medfilt

from .config import HAMPEL_K, HAMPEL_WINDOW


def hampel(x):
    """Return the filtered spectrum and the number of points replaced."""
    med = medfilt(x, HAMPEL_WINDOW)
    res = x - med
    sigma = 1.4826 * np.median(np.abs(res)) + 1e-9
    spikes = np.abs(res) > HAMPEL_K * sigma
    y = x.copy()
    y[spikes] = med[spikes]
    return y, int(spikes.sum())


def noise_sigma(wl, x, lo=400, hi=700):
    """Per-pixel noise (percentage points): robust standard deviation of first differences / sqrt(2)."""
    v = (wl >= lo) & (wl <= hi)
    return float(1.4826 * np.median(np.abs(np.diff(x[v]))) / np.sqrt(2))
