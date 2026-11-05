"""Step 6 - resampling to a common grid: 380-780 nm in 1 nm steps, by linear interpolation."""
import numpy as np

from .config import WL_OUT

GRID = np.arange(WL_OUT[0], WL_OUT[1] + 1, 1.0)


def resample(wl, x):
    return np.interp(GRID, wl, x)
