"""Parameters and species table shared by every step of the Morpho angular reflectance pipeline."""
import numpy as np

# --- wavelength ranges -------------------------------------------------------
WL_OUT = (380, 780)        # output grid of the cleaned spectra, 1 nm steps
WL_WORK = (350, 820)       # range processed before resampling (margin for the filters)
WL_VIS = (400, 700)        # range of the Montgomerie descriptors and of the noise estimate

# --- cleaning ----------------------------------------------------------------
HOT_PIXEL_FRACTION = 0.30  # a pixel is "hot" if it spikes in more than 30 % of the spectra of a species
HOT_PIXEL_MEDIAN = 7       # median-filter window (pixels) used to find spikes
HOT_PIXEL_Z = 5.0          # spike threshold, in robust standard deviations
HAMPEL_WINDOW = 7          # Hampel filter window (pixels)
HAMPEL_K = 6.0             # Hampel threshold, in MAD-based standard deviations
SG_WINDOW = 25             # Savitzky-Golay window (pixels, about 14 nm)
SG_ORDER = 3               # Savitzky-Golay polynomial order

# --- measurement checks (quality control, not colour descriptors) -------------
ROW_DRIFT_MAX_PP = 3.0     # |row offset| above this flags row drift (percentage points)
NEGATIVE_BASELINE_PCT = -2.0
SUSPECT_FACTOR = 5.0       # isolated cell more than 5x brighter than its neighbours

# --- colour ----------------------------------------------------------------
ILLUMINANT = 'D65'
OBSERVER = 'CIE 1931 2 Degree Standard Observer'
EXPOSURE_PERCENTILE = 99   # luminance percentile (all species together) shown as white

# --- goniometer grid --------------------------------------------------------
I_GRID = np.arange(-80, 81, 10)   # incidence angles (deg)
C_GRID = np.arange(-90, 91, 5)    # collection angles (deg)


def in_blind_zone(I, C):
    """Geometries that cannot be measured because the detector blocks the source."""
    return abs(C - I) <= 5 or (abs(I) == 80 and abs(C) == 85 and np.sign(C) == np.sign(I))


# --- species ----------------------------------------------------------------
# (species_id, species name, regular expression matching the raw file names)
# Raw file names keep the original spellings (amphitrion, cisseus, sulkwoskyi); the ids and names use the accepted ones.
SPECIES = [
    ('absoloni', 'Morpho absoloni', 'M. absoloni'),
    ('achilles', 'Morpho achilles', 'M. achilles'),
    ('aega', 'Morpho aega', 'M. aega'),
    ('amathonte', 'Morpho amathonte', 'M. amathonte'),
    ('amphitryon', 'Morpho amphitryon', 'M. amphitrion'),
    ('anaxibia', 'Morpho anaxibia', 'M. anaxibia'),
    ('aurora', 'Morpho aurora', 'M. aurora'),
    ('cisseis', 'Morpho cisseis', 'M. cisseus'),
    ('cypris', 'Morpho cypris', 'M. cypris'),
    ('deidamia', 'Morpho deidamia', 'M. deidamia'),
    ('epistrophus', 'Morpho epistrophus', 'M. epistrophus'),
    ('eugenia', 'Morpho eugenia', 'M. eugenia'),
    ('godarti', 'Morpho godarti', 'M. godarti'),
    ('granadensis', 'Morpho granadensis', 'M. granadensis'),
    ('helenor_p', 'Morpho helenor (p)', 'M. helenor p'),
    ('helenor_unal', 'Morpho helenor (UNAL)', 'M. helenor unal'),
    ('hercules', 'Morpho hercules', 'M. hercules'),
    ('iphitus', 'Morpho iphitus', 'M. iphitus'),
    ('lympharis', 'Morpho lympharis', 'M. lympharis'),
    ('marcus', 'Morpho marcus', 'M. marcus'),
    ('menelaus', 'Morpho menelaus', 'M. menelaus'),
    ('niepelti', 'Morpho niepelti', 'M. niepelti'),
    ('polyphemus', 'Morpho polyphemus', 'M. polyphemus'),
    ('portis', 'Morpho portis', 'M. portis'),
    ('rhetenor', 'Morpho rhetenor', 'M. rhetenor'),
    ('rhodopteron', 'Morpho rhodopteron', 'M. rhodopteron'),
    ('sulkowskyi', 'Morpho sulkowskyi', 'M. sulkwoskyi'),
    ('telemachus', 'Morpho telemachus', 'M. telemachus'),
    # The Spanish sign words are literal source-filename variants and are kept for input matching.
    ('theseus', 'Morpho theseus', 'M. theseus negativo|M. theseus positivo'),
    ('theseus_w', 'Morpho theseus (w)', 'M. theseus w'),
    ('zephyritis', 'Morpho zephyritis', 'M. zephyritis'),
]

# Species measured twice: the remeasurement is the main dataset, the earlier sweep is kept for comparison.
REMEASURED = {'deidamia', 'eugenia', 'lympharis', 'menelaus', 'rhetenor'}

# dataset variants
SINGLE, REMEASUREMENT, PREVIOUS = 'single', 'remeasurement', 'previous_sweep'
