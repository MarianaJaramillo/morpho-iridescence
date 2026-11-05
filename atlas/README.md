# Morpho Angular Atlas: data, scripts and viewer

Angular reflectance of 31 *Morpho* specimens measured with a goniometer and an Avantes spectrometer. This directory contains the data shown in the Morpho Angular Atlas, the processing scripts, and the viewer. The original acquisition workbooks used for the initial cleaning step are not included; the included `raw_spectra/` files are angle-label-corrected data products, not untouched instrument exports.

```
data/        cleaned data and derived tables (generated column names are English; original headers are retained for provenance)
scripts/     cleaning, colour maps, pavo check, Gruson overlay, viewer build
viewer/      the Atlas (single offline file)
```

## 1. Data (`data/`)

### Conventions

| Item | Value |
|---|---|
| Incidence angle `I_deg` | −80° to 80° in 10° steps (angle of the light source from the wing normal) |
| Collection angle `C_deg` | −90° to 90° in 5° steps (angle of the detector from the wing normal) |
| Ideal specular | C = −I |
| Blind zone | \|C − I\| ≤ 5°: the detector blocks the source; never measured |
| Reflectance | % of the white reference. Can exceed 100 % near specular |
| Wavelength | cleaned: 380–780 nm in 1 nm steps; raw: detector pixels in 380–780 nm (~0.57 nm) |
| Acquisition | 500 ms integration, 1 average, smoothing 1 (all measurements) |
| `variant` | `single`: the only sweep of a species · `remeasurement`: main sweep of deidamia, eugenia, lympharis, menelaus and rhetenor · `previous_sweep`: their earlier sweep, kept for comparison |

Species ids use the accepted names (*amphitryon*, *cisseis*, *godarti*, *sulkowskyi*); original acquisition filenames and headers keep their source spellings. These source strings are preserved for provenance and compatibility and must not be renamed without updating and validating the configuration, cleaning rules, indexes, tests, and viewer. The Gruson helenor results are linked only to *Morpho helenor* (UNAL).

### Files

| File | Content |
|---|---|
| `cleaned_spectra/<species_id>.csv` | One row per geometry: `I_deg`, `C_deg`, then `380` … `780` (R in %, 1 nm). Previous sweeps: `<species_id>__previous_sweep.csv` |
| `raw_spectra/<species_id>.csv` | Same layout with the raw detector wavelengths (`380.47`, `381.06`, …). Angles repaired, values untouched |
| `spectra_index.csv` | One row per spectrum: geometry, source file, source column, original header, integration time, averages, spikes removed, noise |
| `metrics_by_geometry.csv` | **Main table.** One row per (species, I, C), 31 species: descriptors, map colours, measurement checks |
| `metrics_previous_sweeps.csv` | Same for the 5 previous sweeps |
| `pavo_descriptors.csv` | B1, B2, B3, H1, S6, S8 computed with the R package pavo |
| `pavo_check.csv` | Largest difference between the Python and pavo descriptors |
| `dataset_summary.csv` | One row per dataset: spectra, missing geometries, hot pixels, noise, row drift, specular peak, tilt |
| `coverage_grid.csv` | Every I × C position per dataset: `measured`, `blind_zone` or `missing` |
| `hot_pixels.csv` | Hot pixels found in each dataset (nm) |
| `cleaning_log.csv` | Every correction: replaced files, header fixes, I/C swaps, copies, duplicates, hot pixels |
| `measurement_checks_log.csv` | Row drift, suspect labels and missing geometries |
| `parameters.json` | Every parameter used, including the colour exposure |
| `gruson/input/` | The Gruson et al. (2019) fit tables (summary, sweeps, points) |
| `gruson/gruson_fits.json` | The same, as read by the viewer |
| `gruson/gruson_alignment.csv` | Shift d used to place the Gruson points on each map, and the correlation with B2 |

### Columns of `metrics_by_geometry.csv`

| Column | Meaning |
|---|---|
| `species_id`, `species`, `variant`, `I_deg`, `C_deg` | Identification |
| `B1` | Total brightness: Σ R over 400–700 nm (%·nm) |
| `B2` | Mean brightness: mean R over 400–700 nm (%) |
| `B3` | Intensity: max R over 400–700 nm (%) |
| `H1` | Hue: wavelength of max R (nm) |
| `S6` | Contrast: max R − min R (percentage points) |
| `S8` | Chroma: (max R − min R) / B2 |
| `perceived_colour_hex` | Colour seen by the eye (D65, 2° observer), sRGB, common exposure for all species |
| `cie_hue_hex` | The same colour with brightness removed (brightest channel = 1) |
| `source_file`, `source_column`, `original_header` | Raw file, column (1 = first data column) and header as typed |
| `integration_time_ms`, `n_averages` | Acquisition settings |
| `n_spikes_removed` | Points replaced by the Hampel filter |
| `noise_sigma_pp` | Per-pixel noise of the spectrum over 400–700 nm (percentage points) |
| `row_offset_pp` | Median B2 of the incidence row away from specular (\|C + I\| ≥ 40°) minus the neighbouring rows |
| `flag_row_drift` | \|`row_offset_pp`\| > 3 pp |
| `flag_negative_baseline` | Median of the spectrum < −2 % (dark-reference offset) |
| `flag_above_white_std` | B3 > 100 % (expected near specular) |
| `flag_suspect_label` | Isolated cell more than 5× brighter than its neighbours, away from specular: probably labelled with the wrong geometry |

B1–S8 are the Montgomerie (2006) variables with the definitions of `pavo::summary.rspec`, computed over 400–700 nm with negative values set to zero (`procspec(fixneg = "zero")`). H1, B3, S6 and S8 are sensitive to noise: in dark, flat cells H1 falls on noise.

## 2. Scripts (`scripts/`)

Python ≥ 3.10 (`python -m pip install -r scripts/requirements.txt`) and, for step 3, R ≥ 4.1 with `pavo` ≥ 2.9. The Python requirements list is not a lockfile; record the resolved package versions when producing or validating a release.

| Step | Script | Does |
|---|---|---|
| 1 | `01_clean_spectra.py` | Cleaning, from the raw xlsx to `cleaned_spectra/`, `raw_spectra/`, `spectra_index.csv`, `cleaning_log.csv` |
| 2 | `02_colour_maps.py` | Descriptors B1–S8, CIE colour (perceived colour, CIE hue), measurement checks |
| 3 | `03_pavo_descriptors.R` | B1, B2, B3, H1, S6, S8 with pavo |
| 4 | `04_check_pavo.py` | Compares steps 2 and 3 |
| 5 | `05_gruson_fits.py` | Reads the Gruson fit tables |
| 6 | `06_gruson_alignment.py` | Places the Gruson points on the maps |
| 7 | `07_build_viewer.py` | Builds the Atlas viewer |

The cleaning and colour code is in `scripts/morpho/`, one module per step:

| Module | Step | Method |
|---|---|---|
| `io_raw.py` | reading | Avantes export saved as xlsx: row 0 names, row 2 integration time, row 3 averages, row 6+ wavelength and R. Corrected files replace the file of the same species and sign |
| `headers.py` | 1. header repair | '+' inside numbers, '+C', letter O for zero, missing space, missing C, I/C swapped in a file, \|C\| > 90 (last two digits), off-grid C assigned to the free slot with the most similar text |
| `dedup.py` | 2. deduplication | Exact copies (hash of values rounded to 4 decimals) dropped; repeated geometry: keep the file whose sign matches I. Final check that every spectrum is still the one of its column |
| `hot_pixels.py` | 3. hot pixels | Residual from a 7-pixel median > 5 robust SD in > 30 % of the spectra of the dataset. Pixel and its two neighbours interpolated |
| `hampel.py` | 4. Hampel filter | Window 7 pixels, threshold 6 robust SD; spikes replaced by the median |
| `smoothing.py` | 5. Savitzky–Golay | Window 25 pixels (~14 nm), polynomial order 3 |
| `resample.py` | 6. resampling | Linear interpolation to 380–780 nm, 1 nm |
| `descriptors.py` | colour descriptors | B1, B2, B3, H1, S6, S8 (Montgomerie 2006 / pavo) |
| `cie_colour.py` | colour maps | CIE XYZ with D65 and the 2° observer (Y = 1 for the white reference) → linear sRGB → perceived colour (÷ 99th percentile of Y of all datasets) and CIE hue (÷ brightest channel) |
| `qa.py` | measurement checks | Row drift, negative baseline, above white, suspect label, tilt |
| `config.py` | | All parameters and the species table |

### Reproduce the processing steps

```bash
cd scripts
# step 1 needs the original xlsx files
python 01_clean_spectra.py  <Spectra folder> <remeasurements folder> <corrections folder> ../data
python 02_colour_maps.py    ../data
Rscript 03_pavo_descriptors.R ../data
python 04_check_pavo.py     ../data
python 05_gruson_fits.py    ../data/gruson/input ../data
python 06_gruson_alignment.py ../data
python 07_build_viewer.py   ../data ../viewer
```

Step 1 requires the original spectrometer XLSX files, which are not included in this repository. Therefore, a clean-room rerun of step 1 is not possible from this distribution. `data/raw_spectra/` holds the spectral values with angle labels repaired; the values are otherwise unprocessed at that stage. Steps 2–7 can run from the included data products. The exact upstream files expected by step 1 are:

* `<Spectra folder>`: the 64 XLSX files from `Spectra.zip`.
* `<remeasurements folder>`: the 10 remeasured XLSX files for deidamia, eugenia, lympharis, menelaus, and rhetenor, described in the source materials as `Datos_angulares_UdeaA`.
* `<corrections folder>`: `M. cisseus negativo paralelo (archivo corregido).xlsx` and `M. amphitrion negativo paralelo (archivo corregido).xlsx`.

The upstream archives or stable download locations, versions, checksums, and redistribution permissions have not been recorded here. Obtain these details from the data owners before claiming that step 1 is reproducible or redistributing those workbooks. Filenames above are reproduced exactly as source identifiers; the Spanish words in them are not translated because the cleaning script matches those names.

## 3. Viewer (`viewer/`)

* `Morpho_Angular_Atlas_offline.html`: double-click to open it in a modern browser. Data, code, and viewer styles are embedded in the file; no server or network connection is required. System font fallbacks are used when IBM Plex fonts are not installed locally.

The viewer build script deletes and recreates `<viewer output folder>/web/`. Do not use a directory containing files that must be preserved as the output folder.

The page has four parts: the explorer (map modes Perceived colour, CIE hue, Brightness B2, Hue H1, Chroma S8; click a cell for its raw and cleaned spectrum), the Gruson model fits per species over the map, the definitions of every quantity, and the cleaning pipeline.

## References

* Montgomerie R. (2006) Analyzing colors. In Hill G.E. & McGraw K.J. (eds.) *Bird Coloration*, vol. 1, Harvard University Press.
* Maia R., Eliason C.M., Bitton P.-P., Doucet S.M., Shawkey M.D. (2013) pavo: an R package for the analysis, visualization and organization of spectral data. *Methods in Ecology and Evolution* 4: 906–913. Maia R., Gruson H., Endler J.A., White T.E. (2019) pavo 2. *Methods in Ecology and Evolution* 10: 1097–1107.
* CIE 15:2018 Colorimetry, 4th edition.
* Gruson H. et al. (2019) Quantitative characterization of iridescent colours in biological studies: a novel method using optical theory. *Interface Focus* 9: 20180049.
* Savitzky A., Golay M.J.E. (1964) *Analytical Chemistry* 36: 1627–1639. Hampel F.R. (1974) *Journal of the American Statistical Association* 69: 383–393.
