"""Colour as seen by the eye: CIE D65 illuminant and CIE 1931 2-degree standard observer (CIE 15:2018).

1. Tristimulus values on the 380-780 nm, 1 nm grid, normalised so that the white reference has Y = 1:
       X = sum R(l) S(l) x(l) / sum S(l) y(l)      (same for Y and Z with y and z), R as a fraction (0-1)
   S = D65 spectral power, x y z = colour-matching functions of the 2-degree observer. Negative R set to zero.
2. Linear sRGB (D65 white point, no chromatic adaptation needed) with negative channels clipped to zero.
3. Perceived colour: linear sRGB divided by one common exposure for all species (the 99th percentile of Y over all
   datasets, previous sweeps included), clipped to 0-1 and gamma-encoded. Brightness is comparable between species and cells.
4. CIE hue: linear sRGB divided by its largest channel (brightest channel = 1), then gamma-encoded. Only the
   chromaticity is left; dark cells become visible, and in very dark cells this colour is mostly noise.
"""
import colour
import numpy as np

from .config import ILLUMINANT, OBSERVER, WL_OUT

_shape = colour.SpectralShape(WL_OUT[0], WL_OUT[1], 1)
_cmfs = colour.MSDS_CMFS[OBSERVER].copy().align(_shape)
_ill = colour.SDS_ILLUMINANTS[ILLUMINANT].copy().align(_shape)
WHITE_POINT = colour.CCS_ILLUMINANTS[OBSERVER][ILLUMINANT]
_K = _ill.values[:, None] * _cmfs.values
_K = _K / (_ill.values * _cmfs.values[:, 1]).sum()


def xyz(spectra_pct):
    """spectra_pct: (n, 401) on the 380-780 nm grid, in %."""
    return np.clip(spectra_pct, 0, None) @ _K / 100.0


def linear_srgb(XYZ):
    return np.clip(colour.XYZ_to_RGB(XYZ, colour.RGB_COLOURSPACES['sRGB'], WHITE_POINT, apply_cctf_encoding=False), 0, None)


def to_hex(rgb01):
    enc = colour.models.eotf_inverse_sRGB(np.clip(rgb01, 0, 1))
    return ['%02x%02x%02x' % tuple((c * 255).round().astype(int)) for c in np.clip(enc, 0, 1)]


def perceived_colour(rgb_lin, exposure_Y):
    return to_hex(np.clip(rgb_lin / exposure_Y, 0, 1))


def cie_hue(rgb_lin):
    mx = rgb_lin.max(1, keepdims=True)
    mx[mx == 0] = 1
    return to_hex(rgb_lin / mx)
