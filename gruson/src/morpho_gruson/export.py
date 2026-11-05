"""Write the tidy, documented output tables (see docs/data_dictionary.md)."""
from pathlib import Path
import numpy as np
import pandas as pd
from .models import brightness_model, hue_model
from .pipeline import _edge_spectrum


def _sweep_label(s):
    return f"{s['session']} #{s['replicate']}" if s["replicate"] else s["session"]


def export_all(R, out_dir):
    out = Path(out_dir); (out / "spectra").mkdir(parents=True, exist_ok=True)
    cfg, final, sessions, pool = R["cfg"], R["final"], R["sessions"], R["pool"]
    sw_b, sw_h, pt_b, pt_h, index = [], [], [], [], []
    wide = {}
    w = next(iter(R["grids"].values()))
    for s in final:
        sp, kind, edge = s["species"], s["kind"], s["method"] == "edge"
        wide.setdefault(sp, {"wavelength_nm": w})
        base = dict(species=sp, session=s["session"], replicate=s["replicate"], n_points=s["n"], R2=round(s["R2"], 5))
        if kind == "brightness":
            sw_b.append(dict(base, method="brightness_mean", Bmax=s["Bmax"], t_deg=s["t"], gammaB_deg=s["gammaB"],
                             constant_span_I_plus_C=s["key"]))
        else:
            sw_h.append(dict(base, method=s["method"], Hmax_nm=s["Hmax"], gammaH=s["gammaH"], constant_I_minus_C=s["key"],
                             note="Hmax = edge position at span 0 (not comparable across species)" if edge else ""))
        sid = f"{sp}|{_sweep_label(s)}|{kind}|key{s['key']}"
        for x, y, I, C, h in zip(s["x"], s["y"], s["I"], s["C"], s["columns"]):
            fit = brightness_model(x, s["Bmax"], s["t"], s["gammaB"]) if kind == "brightness" else hue_model(x, s["Hmax"], s["gammaH"])
            (pt_b if kind == "brightness" else pt_h).append(dict(
                species=sp, session=s["session"], replicate=s["replicate"], I_deg=I, C_deg=C, x_deg=x,
                y_observed=y, y_fitted=fit, residual=y - fit, raw_column=h))
            spec = _edge_spectrum(cfg, sessions[s["session"]][0], sessions[s["session"]][1][h], I, C, pool)[1] if edge \
                else R["spectra"][(s["session"], h)]
            cid = f"{sid}|I{I}C{C}"
            if cid in wide[sp]:
                cid += "|rep2"
            wide[sp][cid] = spec
            index.append(dict(species=sp, column_id=cid, session=s["session"], replicate=s["replicate"], sweep_kind=kind,
                              I_deg=I, C_deg=C, span_I_plus_C=I + C, bisector_deg=(I - C) / 2, raw_column=h,
                              sweep_R2=round(s["R2"], 5)))
    fmt = dict(index=False, encoding="utf-8-sig", float_format="%.6g")
    pd.DataFrame(sw_b).to_csv(out / "sweeps_brightness.csv", **fmt)
    pd.DataFrame(sw_h).to_csv(out / "sweeps_hue.csv", **fmt)
    pd.DataFrame(pt_b).to_csv(out / "fit_points_brightness.csv", **fmt)
    pd.DataFrame(pt_h).to_csv(out / "fit_points_hue.csv", **fmt)
    pd.DataFrame(index).to_csv(out / "spectra" / "INDEX_spectra_columns.csv", **fmt)
    for sp, d in wide.items():
        pd.DataFrame(d).to_csv(out / "spectra" / f"{sp}_spectra.csv", **fmt)
    rows = []
    for sp, r in R["summary"].items():
        e = r["hue_method"] == "edge"
        rows.append(dict(species=sp, hue_method=r["hue_method"], n_brightness=r["Bmax"][2],
                         Bmax_mean=r["Bmax"][0], Bmax_SD=r["Bmax"][1], t_deg_mean=r["t"][0], t_deg_SD=r["t"][1],
                         gammaB_deg_mean=r["gammaB"][0], gammaB_deg_SD=r["gammaB"][1], n_hue=r["gammaH"][2],
                         Hmax_nm_mean=None if e else r["Hmax"][0], Hmax_nm_SD=None if e else r["Hmax"][1],
                         gammaH_mean=r["gammaH"][0], gammaH_SD=r["gammaH"][1]))
    pd.DataFrame(rows).to_csv(out / "summary_by_species.csv", **fmt)
    return len(sw_b), len(sw_h)
