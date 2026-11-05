"""End-to-end Gruson pipeline: raw spreadsheets -> accepted sweeps -> per-species parameters.

Run ``python -m morpho_gruson.cli run`` (see README). Every threshold comes from
``config/pipeline.yaml``.
"""
from collections import defaultdict
from pathlib import Path
import numpy as np
import yaml
from scipy.optimize import curve_fit

from .io import load_sheet, parse_columns, REFPAT, REFPAT2
from .models import brightness_model, hue_model, r_squared
from .signal import smooth, mask_band, nearest_reference, hue_peak_centroid, hue_edge


# ----------------------------------------------------------------------------- loading
def load_config(path):
    return yaml.safe_load(open(path, encoding="utf-8"))


class RawStore:
    """Lazy cache of the raw spreadsheets."""
    def __init__(self, raw_dir):
        self.dir, self._c = Path(raw_dir), {}

    def get(self, fname):
        if fname not in self._c:
            self._c[fname] = load_sheet(self.dir / fname)
        return self._c[fname]


def _checksum(y):
    return round(float(np.nansum(y[300:340])), 2)


def unique_columns(cfg, store):
    """{file: [column names kept]} after removing columns repeated in earlier files.

    Lists, not sets: a set of strings iterates in a different order on every run
    (hash randomisation), which reordered points and spectra in the output tables.
    """
    seen, keep = set(), {}
    for fname in cfg["dedup_order"]:
        _, data = store.get(fname)
        recs, _ = parse_columns(list(data))
        keep[fname] = []
        for pts in recs.values():
            for _, _, _, h in pts:
                k = _checksum(data[h])
                if k in seen:
                    continue
                seen.add(k)
                keep[fname].append(h)
    return keep


def _rename(data, rules):
    out = {}
    for h, v in data.items():
        for r in rules or []:
            if r["mode"] == "substring":
                h = h.replace(r["find"], r["replace"])
            elif r["mode"] == "prefix" and h.startswith(r["find"]):
                h = h.replace(r["find"], r["replace"], 1)
        out[h] = v
    return out


def load_sessions(cfg, store):
    """{session id: (wavelength, {column: raw intensity})} after composition, de-duplication and renames."""
    uniq = unique_columns(cfg, store)
    out = {}
    for s in cfg["sessions"]:
        if "base" in s:
            wl, data = store.get(s["base"]); data = dict(data)
            _, d2 = store.get(s["add_unique"])
            for h in uniq[s["add_unique"]]:
                data[h] = d2[h]
        elif s.get("unique_only"):
            wl, d = store.get(s["file"]); data = {h: d[h] for h in uniq[s["file"]]}
        else:
            wl, data = store.get(s["file"]); data = dict(data)
        out[s["id"]] = (wl, _rename(data, s.get("renames")))
    return out


def reference_pool(cfg, store):
    lo, hi = cfg["signal"]["wavelength_range_nm"]; sg = cfg["signal"]["savgol"]
    pool = {}
    for fname in cfg["reference_pool"]["files"]:
        wl, data = store.get(fname)
        m = (wl >= lo) & (wl <= hi)
        for h in data:
            mm = REFPAT.match(h) or REFPAT2.match(h)
            if not mm:
                continue
            key = (int(mm.group(1)), int(mm.group(2)))
            sm = smooth(data[h], sg["window"], sg["polyorder"])[m]
            pool[key] = (pool[key] + sm) / 2.0 if key in pool else sm
    return pool


# ----------------------------------------------------------------------------- sweeps
def _cluster(items, key):
    """Group items whose integer `key` values differ by <= 1 (a sweep keeps I+C or I-C constant)."""
    ks = sorted({round(key(p)) for p in items}); rep, cur = {}, None
    for k in ks:
        if cur is None or k - cur > 1:
            cur = k
        rep[k] = cur
    out = defaultdict(list)
    for p in items:
        out[rep[round(key(p))]].append(p)
    return out


def normalised_spectra(cfg, sessions, pool):
    """{(session, column): reflectance on the 400-700 nm grid}."""
    lo, hi = cfg["signal"]["wavelength_range_nm"]; sg = cfg["signal"]["savgol"]; band = cfg["signal"]["masked_band_nm"]
    geoms = sorted(pool); spectra, grids = {}, {}
    for sid, (wl, data) in sessions.items():
        m = (wl >= lo) & (wl <= hi); w = wl[m]; grids[sid] = w
        recs, _ = parse_columns(list(data))
        for pts in recs.values():
            for kind, I, C, h in pts:
                sm = smooth(data[h], sg["window"], sg["polyorder"])[m]
                ref = pool[nearest_reference(I, C, geoms)]
                spectra[(sid, h)] = mask_band(w, sm / np.where(ref == 0, np.nan, ref), band)
    return spectra, grids


def build_sweeps(cfg, sessions, spectra, grids):
    """One record per sweep (a set of >= 1 spectra with a constant I+C [brightness] or I-C [hue])."""
    frac = cfg["signal"]["hue_peak_fraction"]; sweeps = []
    for sid, (wl, data) in sessions.items():
        w = grids[sid]
        recs, _ = parse_columns(list(data))
        for (sp, rep), pts in recs.items():
            bp, hp = [], []
            for kind, I, C, h in pts:
                refl = spectra[(sid, h)]
                if kind == "brightness":
                    bp.append(((I - C) / 2.0, float(np.nanmean(refl)), I, C, h))
                else:
                    hp.append((float(I + C), float(hue_peak_centroid(w, refl, frac)), I, C, h))
            for kind, items, keyf in (("brightness", bp, lambda p: p[2] + p[3]), ("hue", hp, lambda p: p[2] - p[3])):
                for k, grp in _cluster(items, keyf).items():
                    sweeps.append(dict(session=sid, species=sp, replicate=rep, kind=kind, key=k, method="peak",
                                       x=[p[0] for p in grp], y=[p[1] for p in grp],
                                       I=[p[2] for p in grp], C=[p[3] for p in grp], columns=[p[4] for p in grp]))
    return sweeps


# ----------------------------------------------------------------------------- fitting
def fit_sweep(sw):
    x, y = np.asarray(sw["x"], float), np.asarray(sw["y"], float)
    sw.update(n=len(x), R2=np.nan)
    if len(x) < 5:
        return sw
    try:
        if sw["kind"] == "brightness":
            po, _ = curve_fit(brightness_model, x, y, p0=[y.max(), x[np.argmax(y)], 10],
                              bounds=([0, -45, 0.5], [np.inf, 45, 60]), maxfev=20000)
            sw.update(Bmax=po[0], t=po[1], gammaB=po[2], R2=r_squared(y, brightness_model(x, *po)))
        else:
            p0 = [y.max(), 0.3 if sw["method"] == "edge" else 0.7]
            po, _ = curve_fit(hue_model, x, y, p0=p0, maxfev=20000)
            sw.update(Hmax=po[0], gammaH=po[1], R2=r_squared(y, hue_model(x, *po)))
    except Exception:
        pass
    return sw


def accepted(sw, acc):
    if sw["n"] < acc["min_points"] or not np.isfinite(sw["R2"]) or sw["R2"] < acc["min_r2"]:
        return False
    if sw["kind"] == "brightness":
        a = acc["brightness"]; lo, hi = a["gammaB_range_deg"]
        return lo < sw["gammaB"] < hi and sw["Bmax"] > a["Bmax_min_exclusive"] and abs(sw["t"]) < a["abs_t_max_deg"]
    return sw["gammaH"] > acc["hue"]["gammaH_min_exclusive"]


def _matches(sw, rule):
    return (sw["species"], sw["kind"], sw["session"], sw["replicate"]) == \
           (rule["species"], rule["kind"], rule["session"], str(rule["sweep"]))


def apply_exclusions(cfg, sweeps):
    for rule in cfg.get("point_exclusions", []):
        for sw in sweeps:
            if _matches(sw, rule) and len(sw["x"]) >= 5:
                keep = [i for i, x in enumerate(sw["x"]) if abs(x - rule["bisector_deg"]) > 1e-6]
                for f in ("x", "y", "I", "C", "columns"):
                    sw[f] = [sw[f][i] for i in keep]
    for sw in sweeps:
        fit_sweep(sw)
    for sw in sweeps:
        sw["manual_exclusion"] = any(_matches(sw, r) for r in cfg["manual_exclusions"])
        sw["accepted"] = accepted(sw, cfg["acceptance"]) and not sw["manual_exclusion"]
    return sweeps


def _edge_spectrum(cfg, wl, raw, I, C, pool):
    """Normalised spectrum used by the edge method.

    Legacy detail kept so published numbers are reproduced exactly: when two white references are equally
    close in (I, C), the edge method takes the one that appears first in the pool (insertion order), whereas
    the peak/brightness path takes the lowest (I, C). Affects only ties, changes gammaH by <1e-3.
    """
    lo, hi = cfg["signal"]["wavelength_range_nm"]; sg = cfg["signal"]["savgol"]
    m = (wl >= lo) & (wl <= hi); w = wl[m]
    g = min(pool, key=lambda q: abs(I - q[0]) + abs(C - q[1]))
    sm = smooth(raw, sg["window"], sg["polyorder"])[m]
    return w, mask_band(w, sm / np.where(pool[g] == 0, np.nan, pool[g]), cfg["signal"]["masked_band_nm"])


def edge_sweeps(cfg, sessions, spectra, grids, pool):
    """Hue sweeps for species without a spectral peak: position of the spectral edge."""
    ec = cfg["edge_method"]; out = []
    for sid, (wl, data) in sessions.items():
        wl_s = sessions[sid][0]; recs, _ = parse_columns(list(data))
        for (sp, rep), pts in recs.items():
            if sp not in ec["species"]:
                continue
            hp = [(I, C, h) for kind, I, C, h in pts if kind == "hue"]
            for k, grp in _cluster(hp, lambda p: p[0] - p[1]).items():
                xs, ys, cols = [], [], []
                for I, C, h in grp:
                    w, spec = _edge_spectrum(cfg, wl_s, data[h], I, C, pool)
                    e = hue_edge(w, spec, ec["plateau_nm"], ec["floor_nm"])
                    if np.isfinite(e):
                        xs.append(float(I + C)); ys.append(float(e)); cols.append((I, C, h))
                if len(xs) < ec["min_points"]:
                    continue
                sw = dict(session=sid, species=sp, replicate=rep, kind="hue", key=k, method="edge", x=xs, y=ys,
                          I=[c[0] for c in cols], C=[c[1] for c in cols], columns=[c[2] for c in cols])
                fit_sweep(sw)
                if np.isfinite(sw["R2"]) and sw["R2"] >= ec["min_r2"]:
                    sw.update(accepted=True, manual_exclusion=False); out.append(sw)
    return out


def select_final(cfg, sweeps, edge):
    edge_sp = set(cfg["edge_method"]["species"])
    final = [s for s in sweeps if s["accepted"] and not (s["kind"] == "hue" and s["species"] in edge_sp)] + edge
    order = {s["id"]: i for i, s in enumerate(cfg["sessions"])}
    final.sort(key=lambda s: (s["species"], s["kind"], order[s["session"]], s["replicate"], s["key"]))
    return final


def species_summary(final):
    """Mean, SD (ddof=1) and n of each parameter across accepted sweeps."""
    def ms(v):
        v = np.asarray(v, float)
        return (v.mean(), v.std(ddof=1) if len(v) > 1 else np.nan, len(v)) if len(v) else (np.nan, np.nan, 0)
    rows = {}
    for sp in sorted({s["species"] for s in final}):
        b = [s for s in final if s["species"] == sp and s["kind"] == "brightness"]
        h = [s for s in final if s["species"] == sp and s["kind"] == "hue"]
        row = {"hue_method": h[0]["method"] if h else ""}
        for p in ("Bmax", "t", "gammaB"):
            row[p] = ms([s[p] for s in b])
        for p in ("Hmax", "gammaH"):
            row[p] = ms([s[p] for s in h])
        rows[sp] = row
    return rows


def run(config_path, raw_dir):
    cfg = load_config(config_path); store = RawStore(raw_dir)
    sessions = load_sessions(cfg, store); pool = reference_pool(cfg, store)
    spectra, grids = normalised_spectra(cfg, sessions, pool)
    sweeps = apply_exclusions(cfg, build_sweeps(cfg, sessions, spectra, grids))
    edge = edge_sweeps(cfg, sessions, spectra, grids, pool)
    final = select_final(cfg, sweeps, edge)
    return dict(cfg=cfg, sessions=sessions, pool=pool, spectra=spectra, grids=grids,
                all_sweeps=sweeps, final=final, summary=species_summary(final))
