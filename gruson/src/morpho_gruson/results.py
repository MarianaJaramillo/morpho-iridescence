"""Helpers to load the processed tables (data/processed/*.csv)."""
from pathlib import Path
import pandas as pd

def _csv(p):
    return pd.read_csv(p, encoding="utf-8-sig", keep_default_na=False, na_values=[""])

def load_processed(processed_dir="data/processed"):
    d = Path(processed_dir)
    out = {n: _csv(d / f"{n}.csv") for n in ("summary_by_species", "sweeps_brightness", "sweeps_hue",
                                             "fit_points_brightness", "fit_points_hue")}
    for k in ("sweeps_brightness", "sweeps_hue", "fit_points_brightness", "fit_points_hue"):
        out[k]["replicate"] = out[k]["replicate"].fillna("").astype(str).str.replace(r"\.0$", "", regex=True)
    return out
