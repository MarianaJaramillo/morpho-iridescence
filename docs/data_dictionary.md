# Data dictionary

All processed tables are UTF-8 with BOM (open directly in Excel), one row per item, snake_case columns.
Angles are in degrees. Wavelengths in nm.

## Raw data: `gruson/data/raw/*.xlsx`
One spreadsheet per measurement session, as exported from the spectrometer software. Column A = wavelength; every other column is one spectrum. File names follow `<YYYY-MM-DD>_<role>.xlsx` (date of the session or export, sorted chronologically). Column headers inside the files are retained exactly as received from the instrument operators, including Spanish words and historical spelling errors; they are source identifiers, not translated project terminology, and changing them can break the pipeline configuration or parser.

| File | Session id | Notes | Original file name |
|---|---|---|---|
| `2026-08-31_S1_main_dataset.xlsx` | S1 | main early dataset + 12 white references | `ultima_version_datos_gruson.xlsx` |
| `2026-08-31_older_export.xlsx` | de-dup only | older export; repeated columns are used once | `31agosto.xlsx` |
| `2026-09-09_older_export.xlsx` | S1 | older export; its unique columns are added to S1 | `9_septiembre.xlsx` |
| `2026-09-02_session.xlsx` | 2026-09-02 | | `2_septiembre_2.xlsx` |
| `2026-09-25_session.xlsx` | 2026-09-25 | | `25_de_septiembre.xlsx` |
| `2026-09-30_session.xlsx` | 2026-09-30 | also provides 27 white references | `30sep.xlsx` |
| `2026-10-07_session.xlsx` | 2026-10-07 | | `7_de_octubre.xlsx` |
| `2026-10-09_session.xlsx` | 2026-10-09 | | `9_de_octubre.xlsx` |
| `2026-10-15_session.xlsx` | 2026-10-15 | | `15_de_octubre.xlsx` |
| `2026-10-19_session.xlsx` | 2026-10-19 | | `19_de_octubre.xlsx` |
| `2026-10-21_session.xlsx` | 2026-10-21 | | `21_octubre.xlsx` |

**Column naming**: `<species><replicate>_<type>_I<I>C<C>.RAW8`, e.g. `aega9_brillo_I13C27.RAW8` = *M. aega*, replicate 9, brightness sweep, I = 13°, C = 27°. Older files use `M. aega Brillo I13 C27.RAW8`. In original headers, `brillo` means brightness and `matiz` means hue; the parser maps these source labels to English internal values. White references: `Ref I20 C20.RAW8` or `ref_I20C20.RAW8`.
Typos corrected in code (`config/pipeline.yaml`): `sulkoeskyi→sulkowskyi`, `aegenia5→eugenia5`, `aega8_brillo_I12C328→I12C28`, `rhodopteron_matiz→rhodopteron5_matiz` (15 Oct file).

## Processed data: `gruson/data/processed/`
Contain the accepted sweeps only (87 brightness, 97 hue).

### `summary_by_species.csv`: one row per species (20)
| Column | Meaning |
|---|---|
| `species` | species epithet (genus *Morpho*) |
| `hue_method` | `peak` or `edge` |
| `n_brightness`, `n_hue` | accepted sweeps |
| `Bmax_mean`, `Bmax_SD` | peak brightness (relative reflectance) |
| `t_deg_mean`, `t_deg_SD` | tilt: bisector of maximal brightness (°) |
| `gammaB_deg_mean`, `gammaB_deg_SD` | angular width of the brightness lobe (°) |
| `Hmax_nm_mean`, `Hmax_nm_SD` | peak wavelength at span 0 (empty for `edge` species) |
| `gammaH_mean`, `gammaH_SD` | hue change with span (dimensionless) |

### `sweeps_brightness.csv` / `sweeps_hue.csv`: one row per accepted sweep
`species`, `session` (ids in `gruson/config/pipeline.yaml`), `replicate` (number in the raw column names; empty if none), `n_points`, `R2`, `method` (`brightness_mean`, `peak`, `edge`), parameters (`Bmax, t_deg, gammaB_deg` | `Hmax_nm, gammaH`), `constant_span_I_plus_C` | `constant_I_minus_C` (the held-constant geometry), `note`.

### `fit_points_brightness.csv` / `fit_points_hue.csv`: one row per spectrum used in a fit
`species, session, replicate, I_deg, C_deg`, `x_deg` (bisector (I−C)/2 for brightness; span I+C for hue), `y_observed` (mean reflectance | peak or edge wavelength in nm), `y_fitted`, `residual`, `raw_column` (original column name).

### `spectra/<species>_spectra.csv`: cleaned spectra
Rows = wavelength (400–700 nm, 520 points, column `wavelength_nm`); columns = spectra with ids
`<species>|<session>[ #<replicate>]|<brightness|hue>|key<held-constant value>|I<I>C<C>` (suffix `|rep2` if the same geometry appears twice in one sweep, as in the *M. godartii* S1 hue sweep, which has 14 points).
`spectra/INDEX_spectra_columns.csv` gives, for every column id: species, session, replicate, sweep kind, I, C, span, bisector, original column name, sweep R².

**Note.** *M. anaxibia* has one accepted brightness sweep (no SD). For *eugenia*, *marcus* and *anaxibia* only γH is comparable (edge method).
