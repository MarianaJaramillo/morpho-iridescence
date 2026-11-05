"""Step 5 - Savitzky-Golay smoothing: window 25 pixels (about 14 nm at the detector's ~0.57 nm pixel), order 3."""
from scipy.signal import savgol_filter

from .config import SG_ORDER, SG_WINDOW


def savitzky_golay(x):
    return savgol_filter(x, SG_WINDOW, SG_ORDER)
