"""Consistency checks for the published Atlas data products (no recomputation of the pipeline)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[1] / "data"
KEY = ["species_id", "variant", "I_deg", "C_deg"]
I_GRID = set(range(-80, 81, 10))
C_GRID = set(range(-90, 91, 5))


@pytest.fixture(scope="module")
def tables():
    read = lambda name: pd.read_csv(DATA / name)
    return {
        "main": read("metrics_by_geometry.csv"),
        "previous": read("metrics_previous_sweeps.csv"),
        "index": read("spectra_index.csv"),
        "coverage": read("coverage_grid.csv"),
        "summary": read("dataset_summary.csv"),
    }


def test_documented_dataset_sizes(tables):
    main, previous = tables["main"], tables["previous"]
    assert main.species_id.nunique() == 32
    assert len(main) == 18483
    assert previous.species_id.nunique() == 5
    assert len(previous) == 2888
    assert set(main.variant) == {"single", "remeasurement"}
    assert set(previous.variant) == {"previous_sweep"}


def test_geometry_keys_are_unique_and_on_the_goniometer_grid(tables):
    for name in ("main", "previous"):
        df = tables[name]
        assert not df.duplicated(KEY).any(), name
        assert set(df.I_deg) <= I_GRID, name
        assert set(df.C_deg) <= C_GRID, name


def test_spectra_index_matches_metrics(tables):
    metrics = pd.concat([tables["main"], tables["previous"]])[KEY]
    index = tables["index"][KEY]
    assert len(index) == len(metrics)
    merged = metrics.merge(index, on=KEY, how="outer", indicator=True)
    assert (merged["_merge"] == "both").all()


def test_dataset_summary_agrees_with_metrics_and_coverage(tables):
    counts = pd.concat([tables["main"], tables["previous"]]).groupby(["species_id", "variant"]).size()
    summary = tables["summary"].set_index(["species_id", "variant"])
    assert len(summary) == len(counts) == 37
    assert (summary.n_spectra == counts.reindex(summary.index)).all()
    coverage = tables["coverage"]
    missing = coverage[coverage.status == "missing"].groupby(["species_id", "variant"]).size()
    assert (summary.n_missing == missing.reindex(summary.index).fillna(0)).all()
    assert set(coverage.status) == {"measured", "blind_zone", "missing"}


def test_measured_cells_in_the_blind_zone_are_the_known_exceptions(tables):
    # Documented as never measured, yet these 9 cells carry spectra; origin to be confirmed by the data owners.
    known = {
        ("aega", "single", -80, -85), ("anaxibia", "single", -10, -15), ("anaxibia", "single", 10, 10),
        ("anaxibia", "single", 40, 35), ("eugenia", "previous_sweep", 10, 5), ("godarti", "single", -20, -20),
        ("godarti", "single", -20, -15), ("lympharis", "previous_sweep", -10, -5), ("zephyritis", "single", 20, 15),
    }
    coverage = tables["coverage"]
    measured = coverage[coverage.status == "measured"]
    inside = measured[(measured.C_deg - measured.I_deg).abs() <= 5]
    assert set(map(tuple, inside[KEY].to_numpy().tolist())) == known


def test_measured_cells_match_metrics(tables):
    measured = tables["coverage"].query("status == 'measured'")
    metrics = pd.concat([tables["main"], tables["previous"]])[KEY]
    merged = metrics.merge(measured[KEY], on=KEY, how="outer", indicator=True)
    assert (merged["_merge"] == "both").all()


def test_descriptors_and_colours_are_well_formed(tables):
    for name in ("main", "previous"):
        df = tables[name]
        for col in ("B1", "B2", "B3", "H1"):
            assert np.isfinite(df[col]).all(), col
        assert df.H1.between(400, 700).all()
        for col in ("perceived_colour_hex", "cie_hue_hex"):
            assert df[col].astype(str).str.fullmatch(r"[0-9a-fA-F]{6}").all(), col


@pytest.mark.parametrize("folder", ["cleaned_spectra", "raw_spectra"])
def test_spectra_files_follow_the_metrics_geometry(tables, folder):
    metrics = pd.concat([tables["main"], tables["previous"]])
    for (species_id, variant), group in metrics.groupby(["species_id", "variant"]):
        suffix = "__previous_sweep" if variant == "previous_sweep" else ""
        spectra = pd.read_csv(DATA / folder / f"{species_id}{suffix}.csv")
        expected = group.sort_values(["I_deg", "C_deg"])
        assert list(spectra.columns[:2]) == ["I_deg", "C_deg"]
        assert (spectra.I_deg.to_numpy() == expected.I_deg.to_numpy()).all(), (folder, species_id, variant)
        assert (spectra.C_deg.to_numpy() == expected.C_deg.to_numpy()).all(), (folder, species_id, variant)


def test_cleaned_spectra_use_the_documented_wavelength_grid():
    wavelengths = [int(float(c)) for c in pd.read_csv(DATA / "cleaned_spectra" / "aega.csv", nrows=1).columns[2:]]
    assert wavelengths == list(range(380, 781))
