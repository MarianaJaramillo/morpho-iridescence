"""Gruson et al. (2019) functional forms and small numerical helpers."""
import numpy as np


def brightness_model(bisector_deg, Bmax, t, gammaB):
    """B(x) = Bmax * exp(-(x - t)^2 / (2 gammaB^2)),  x = bisector (I - C) / 2 in degrees."""
    return Bmax * np.exp(-((bisector_deg - t) ** 2) / (2 * gammaB ** 2))


def hue_model(span_deg, Hmax, gammaH):
    """H(s) = Hmax * cos(gammaH * s),  s = span I + C in degrees (gammaH is dimensionless: per degree)."""
    return Hmax * np.cos(np.deg2rad(gammaH * span_deg))


def r_squared(y, y_fit):
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan


def hyperbola(gammaB, k):
    """Empirical relation Bmax = k / gammaB (area under the Gaussian lobe ~ constant)."""
    return k / gammaB
