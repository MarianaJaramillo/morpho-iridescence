# Methods: the Gruson protocol applied to *Morpho* wings

## 1. Measurement geometry
A goniometer holds an **illumination** fibre at angle **I** and a **collection** fibre at angle **C** (degrees, relative to the sample normal, both in the goniometer's rotation plane, with the scale ridges aligned to that plane). Two derived quantities define the sweeps:

- **bisector** = (I − C) / 2 — the direction the sample "must face" for specular reflection;
- **span** = I + C — the angle between the fibres.

Gruson et al. (2019) showed that, for thin-film-like iridescence, brightness depends mainly on the bisector and hue mainly on the span. Hence two 1-D sweeps are enough:

| Sweep | Held constant | Varied | Quantity extracted |
|---|---|---|---|
| **Brightness** | span I + C | bisector (I − C)/2 | Bmax, t, γB |
| **Hue** | I − C (bisector) | span I + C | Hmax, γH |

(The original instrument exports use the labels `brillo` for brightness and `matiz` for hue. These source labels are intentionally preserved for parsing; output tables use the English values `brightness` and `hue`.)

## 2. Spectral processing (every spectrum)
1. Raw intensity exported by the spectrometer software; analysis window **400–700 nm** (520 points).
2. **Savitzky–Golay** smoothing (window 151 points, polynomial order 3). Missing values are linearly interpolated first.
3. **White-reference normalisation**: divide by the white reference measured at the geometry closest in (I, C) (L1 distance), taken from a pool of **36 references** merged from two sessions (a geometry measured in both is averaged). Values are *relative* reflectance and can be slightly < 0 or > 1.
4. **Artefact mask**: 523–534 nm is replaced by linear interpolation (instrument/lamp artefact).

## 3. Descriptors
- **Brightness** of a spectrum = mean reflectance over 400–700 nm.
- **Hue** of a spectrum = intensity-weighted centroid of the wavelengths within the top 10 % of the peak height (peak method).
- **Edge method** (*M. eugenia*, *M. marcus*, *M. anaxibia*): these spectra have no peak in 400–700 nm (reflectance falls monotonically from the blue plateau). Hue is the wavelength at the half-fall between the plateau (median 420–470 nm) and the floor (median 620–670 nm). Only **γH** is comparable with the peak method; the "Hmax" of these species is an edge position and is not reported.

## 4. Models
- Brightness: **B(x) = Bmax · exp(−(x − t)² / (2 γB²))**, x = bisector. Bmax = peak brightness, t = tilt (bisector of maximum), γB = angular width (°).
- Hue: **H(s) = Hmax · cos(γH · s)**, s = span in degrees; γH is dimensionless (per degree). Larger γH = stronger colour shift with angle.

Fits are non-linear least squares (`scipy.optimize.curve_fit`), brightness bounds t ∈ [−45, 45], γB ∈ [0.5, 60].

## 5. Acceptance of a sweep
A fitted sweep is retained if all of the following hold (`config/pipeline.yaml`); `processed/` contains only retained sweeps:

| Criterion | Value |
|---|---|
| Points in the sweep | ≥ 5 |
| R² | ≥ 0.75 |
| γB | 0.5 < γB < 59 ° |
| Bmax | > 0 |
| \|t\| | < 44 ° |
| γH | > 0.05 |

Edge-method sweeps need ≥ 5 points and R² ≥ 0.75.

## 6. Species level
Species parameters are the **mean ± SD (ddof = 1)** across accepted sweeps (replicates from different sessions are treated as independent sweeps). A species with a single sweep has no SD.

## 7. Additional analyses shown in the presentation
- **Specular vs diffuse**: Bmax against γB per species; hyperbola Bmax = k/γB (k = 7.17, R² = 0.77). It has the best small-sample criteria (AICc, leave-one-out error) among six alternatives we tried, and a physical reading: the lobe area ∝ Bmax·γB is roughly constant.
- **Clustering**: Ward linkage on standardised (ln γB, γH), k = 3 (silhouette 0.525; k = 2: 0.390, k = 4: 0.452).
- **Correlation** γB–γH: Spearman ρ = −0.12 (p = 0.62); no relation detected.

## References
Gruson H., Andraud C., Daney de Marcillac W., Berthier S., Elias M., Gomez D. (2019) Quantitative characterization of iridescent colours in biological studies: a novel method using optical theory. *Interface Focus* 9: 20180049. doi:10.1098/rsfs.2018.0049
