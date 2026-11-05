"""The pipeline must reproduce the species table shown in the presentation (atlas/data/gruson/input/summary_by_species.csv)."""
import warnings
from pathlib import Path
import numpy as np, pandas as pd, pytest
from morpho_gruson.io import parse_columns
from morpho_gruson.pipeline import run

ROOT = Path(__file__).resolve().parents[1]
warnings.filterwarnings("ignore")


@pytest.fixture(scope="module")
def result():
    return run(ROOT / "config/pipeline.yaml", ROOT / "data/raw")


def test_original_sweep_labels_map_to_english():
    raw_headers = ["M. aega8 Brillo I12 C28.RAW8", "aega8_matiz_I5C35.RAW8"]

    records, references = parse_columns(raw_headers)

    assert [record[0] for record in records[("aega", "8")]] == ["brightness", "hue"]
    assert [record[3] for record in records[("aega", "8")]] == raw_headers
    assert references == []


def test_sweep_counts(result):
    f = result["final"]
    assert sum(s["kind"] == "brightness" for s in f) == 87
    assert sum(s["kind"] == "hue" for s in f) == 97
    assert len(result["summary"]) == 20
    assert len(result["pool"]) == 36


def test_species_parameters_match_reference(result):
    ref = pd.read_csv(ROOT.parent / "atlas/data/gruson/input/summary_by_species.csv", encoding="utf-8-sig").set_index("species")
    cols = {"Bmax": "Bmax_mean", "t": "t_deg_mean", "gammaB": "gammaB_deg_mean", "gammaH": "gammaH_mean", "Hmax": "Hmax_nm_mean"}
    for sp, row in result["summary"].items():
        assert row["Bmax"][2] == ref.loc[sp, "n_brightness"] and row["gammaH"][2] == ref.loc[sp, "n_hue"]
        for p, c in cols.items():
            if p == "Hmax" and row["hue_method"] == "edge":
                continue
            assert row[p][0] == pytest.approx(ref.loc[sp, c], rel=1e-4, abs=1e-5), (sp, p)
