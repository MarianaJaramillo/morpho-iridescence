"""Figures of the summary presentation, built only from data/processed/*.csv."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import ConvexHull
from sklearn.metrics import silhouette_score
from adjustText import adjust_text
from .models import brightness_model, hue_model, hyperbola, r_squared

RED, BLU, GRY, ORG, GRN, TEAL = "#C0392B", "#2F5597", "#8A94A3", "#E69F00", "#009E73", "#0072B2"
SPECULAR = ["rhodopteron", "aega", "zephyritis"]; DIFFUSE = ["achilles", "helenor", "deidamia"]
IRIDESCENT = ["menelaus", "aega", "portis", "zephyritis"]; STABLE = ["rhetenor", "deidamia"]
name = lambda s: "M. " + s


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12})

def _clean(a):
    for k in ("top", "right"):
        a.spines[k].set_visible(False)

def _summary_with_sd(R):
    """Species with at least two accepted brightness sweeps (an SD exists), as in the presentation."""
    s = R["summary_by_species"].set_index("species")
    return s[s["gammaB_deg_SD"].notna()]


def sweep_examples(R, species=("aega", "helenor"), cmaps=(plt.cm.Reds, plt.cm.Blues)):
    """All accepted sweeps of two species: brightness (left) and hue (right), with Gaussian / cosine fits."""
    S = R["summary_by_species"].set_index("species")
    fig, ax = plt.subplots(len(species), 2, figsize=(12.5, 3.6 * len(species)))
    for i, (sp, cm) in enumerate(zip(species, cmaps)):
        for j, kind in enumerate(("brightness", "hue")):
            a = ax[i, j]; sw = R[f"sweeps_{kind}"]; sw = sw[sw.species == sp].reset_index(drop=True)
            pts = R[f"fit_points_{kind}"]; n = len(sw)
            for k, r in sw.iterrows():
                c = cm(0.45 + 0.5 * k / max(n - 1, 1))
                p = pts[(pts.species == sp) & (pts.session == r.session) & (pts.replicate == r.replicate)]
                x, y = p.x_deg.values, p.y_observed.values
                if kind == "brightness":
                    xx = np.linspace(x.min() - 3, x.max() + 3, 300); yy = brightness_model(xx, r.Bmax, r.t_deg, r.gammaB_deg)
                else:
                    xx = np.linspace(0, x.max() + 3, 300); yy = hue_model(xx, r.Hmax_nm, r.gammaH)
                a.plot(xx, yy, color=c, lw=1.8, alpha=.9); a.scatter(x, y, color=c, s=34, edgecolor="white", lw=.5, zorder=3)
            s = S.loc[sp]; r2 = sw.R2
            if kind == "brightness":
                a.set_xlabel("Bisector (I−C)/2 (°)"); a.set_ylabel("Brightness (mean reflectance)")
                txt = f"{n} sweeps, $R^2$ = {r2.min():.3f}–{r2.max():.3f}\n$B_{{max}}$ = {s.Bmax_mean:.2f} ± {s.Bmax_SD:.2f}\n$\\gamma_B$ = {s.gammaB_deg_mean:.1f} ± {s.gammaB_deg_SD:.1f}°"
            else:
                a.set_xlabel("Span I+C (°)"); a.set_ylabel("Hue: peak wavelength (nm)")
                txt = f"{n} sweeps, $R^2$ = {r2.min():.3f}–{r2.max():.3f}\n$H_{{max}}$ = {s.Hmax_nm_mean:.0f} ± {s.Hmax_nm_SD:.0f} nm\n$\\gamma_H$ = {s.gammaH_mean:.2f} ± {s.gammaH_SD:.2f}"
            a.text(.97, .95, txt, transform=a.transAxes, ha="right", va="top", fontsize=11, bbox=dict(fc="white", ec="#ccc", boxstyle="round"))
            a.set_title(f"{name(sp)} — {kind}", loc="left", fontweight="bold"); a.grid(alpha=.25); _clean(a)
    fig.tight_layout(); return fig


def specular_vs_diffuse(R):
    """Bmax against gammaB (species means ± SD) with the hyperbola Bmax = k / gammaB."""
    s = _summary_with_sd(R); x, y = s.gammaB_deg_mean.values, s.Bmax_mean.values
    k = curve_fit(hyperbola, x, y, p0=[7])[0][0]; r2 = r_squared(y, hyperbola(x, k)); xs = np.linspace(2.3, 12.8, 200)
    fig, a = plt.subplots(figsize=(10.5, 6.6)); a.plot(xs, k / xs, color="#444", lw=1.6, ls="--", zorder=1); tx = []
    for sp, d in s.iterrows():
        c = RED if sp in SPECULAR else BLU if sp in DIFFUSE else GRY
        a.errorbar(d.gammaB_deg_mean, d.Bmax_mean, xerr=d.gammaB_deg_SD, yerr=d.Bmax_SD, fmt="o", color=c, ms=8, capsize=2, lw=1.1, mec="white")
        tx.append(a.text(d.gammaB_deg_mean, d.Bmax_mean, name(sp), style="italic", fontsize=10.5))
    adjust_text(tx, ax=a, arrowprops=dict(arrowstyle="-", color="grey", lw=.5))
    a.set_xlabel("$\\gamma_B$ — angular width of the reflected lobe (°)   [narrow = directional]"); a.set_ylabel("$B_{max}$ — peak brightness")
    a.text(.97, .95, f"$B_{{max}}$ = {k:.2f} / $\\gamma_B$   ($R^2$ = {r2:.2f})", transform=a.transAxes, ha="right", va="top", fontsize=12, color="#444")
    a.annotate("SPECULAR\n(bright, narrow)", xy=(2.9, 3.15), xytext=(5.2, 3.3), color=RED, fontweight="bold", fontsize=12, arrowprops=dict(arrowstyle="->", color=RED))
    a.annotate("DIFFUSE\n(dim, wide)", xy=(11.4, 0.8), xytext=(9.4, 1.9), color=BLU, fontweight="bold", fontsize=12, arrowprops=dict(arrowstyle="->", color=BLU))
    a.grid(alpha=.25); _clean(a); fig.tight_layout(); return fig, dict(k=k, R2=r2)


def gammaH_bars(R):
    """gammaH per species (mean ± SD); hatched bars = hue from the spectral edge."""
    s = R["summary_by_species"].sort_values("gammaH_mean", ascending=False).reset_index(drop=True)
    fig, a = plt.subplots(figsize=(8.2, 7.4))
    for i, r in s.iterrows():
        ed = r.hue_method == "edge"; c = ORG if r.species in IRIDESCENT else TEAL if r.species in STABLE else GRY
        a.barh(i, r.gammaH_mean, color="white" if ed else c, edgecolor=c, hatch="///" if ed else None, height=.75, lw=1.3)
        a.errorbar(r.gammaH_mean, i, xerr=r.gammaH_SD, color="#222", capsize=2.5, lw=1)
    a.set_yticks(range(len(s))); a.set_yticklabels([name(x) for x in s.species], style="italic", fontsize=11); a.invert_yaxis(); a.set_xlim(.2, .46)
    a.set_xlabel("$\\gamma_H$ — hue change with viewing angle"); a.grid(axis="x", alpha=.25); _clean(a)
    a.legend(handles=[Patch(color=ORG, label="Most iridescent"), Patch(color=TEAL, label="Most stable hue"),
                      Patch(facecolor="white", edgecolor="#555", hatch="///", label="Hue from spectral edge")], frameon=False, fontsize=11, loc="lower right")
    fig.tight_layout(); return fig


def gammaB_vs_gammaH(R, k=3):
    """gammaB vs gammaH with Ward clustering on (ln gammaB, gammaH), standardised."""
    s = _summary_with_sd(R); sp = list(s.index)
    gb, gh = s.gammaB_deg_mean.values, s.gammaH_mean.values
    Z = np.column_stack([np.log(gb), gh]); Z = (Z - Z.mean(0)) / Z.std(0); L = linkage(Z, "ward")
    sil = {kk: silhouette_score(Z, fcluster(L, kk, "maxclust")) for kk in (2, 3, 4)}
    cl = fcluster(L, k, "maxclust"); mg = {c: gb[cl == c].mean() for c in set(cl)}
    order = sorted(mg, key=lambda c: -mg[c]); cl = np.array([order.index(c) for c in cl])
    COL = [ORG, TEAL, GRN]; LAB = ["Wide lobe (diffuse)", "Narrow lobe, stable hue", "Narrow lobe, iridescent"]
    if gh[cl == 1].mean() > gh[cl == 2].mean():
        COL[1], COL[2] = COL[2], COL[1]; LAB[1], LAB[2] = LAB[2], LAB[1]
    fig, a = plt.subplots(figsize=(10, 6.4)); tx = []
    for i in range(k):
        m = cl == i; pts = np.column_stack([gb[m], gh[m]])
        if len(pts) >= 3:
            h = ConvexHull(pts); a.fill(pts[h.vertices, 0], pts[h.vertices, 1], color=COL[i], alpha=.13, lw=0)
    for sp_, c in zip(sp, cl):
        d = s.loc[sp_]; ed = d.hue_method == "edge"
        a.errorbar(d.gammaB_deg_mean, d.gammaH_mean, xerr=d.gammaB_deg_SD, yerr=d.gammaH_SD, fmt="s" if ed else "o",
                   mfc="white" if ed else COL[c], color=COL[c], mec=COL[c], mew=1.5, ms=8, capsize=2, lw=1)
        tx.append(a.text(d.gammaB_deg_mean, d.gammaH_mean, name(sp_), style="italic", fontsize=10))
    adjust_text(tx, ax=a, arrowprops=dict(arrowstyle="-", color="grey", lw=.5))
    a.set_xlabel("$\\gamma_B$ — angular width of brightness (°)"); a.set_ylabel("$\\gamma_H$ — hue change with angle"); a.grid(alpha=.25); _clean(a)
    a.legend(handles=[Patch(color=COL[i], label=LAB[i]) for i in range(k)] + [plt.Line2D([], [], marker="s", mfc="white", color="#555", ls="", label="$\\gamma_H$ from spectral edge")],
             frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    groups = {LAB[i]: [x for x, c in zip(sp, cl) if c == i] for i in range(k)}
    rho = spearmanr(gb, gh)
    return fig, dict(groups=groups, silhouette=sil, spearman_rho=rho.statistic, spearman_p=rho.pvalue)


def results_table(R):
    """Final table: mean ± SD (n) per species, formatted as in the presentation."""
    import pandas as pd
    f = lambda m, s, d: f"{m:.{d}f} ± {s:.{d}f}" if np.isfinite(s) else f"{m:.{d}f}"
    rows = []
    for _, r in R["summary_by_species"].iterrows():
        ed = r.hue_method == "edge"
        rows.append({"Species": name(r.species), "n (B/H)": f"{r.n_brightness}/{r.n_hue}",
                     "Bmax": f(r.Bmax_mean, r.Bmax_SD, 2), "t (°)": f(r.t_deg_mean, r.t_deg_SD, 1), "γB (°)": f(r.gammaB_deg_mean, r.gammaB_deg_SD, 1),
                     "Hmax (nm)": "—" if ed else f(r.Hmax_nm_mean, r.Hmax_nm_SD, 0), "γH": f(r.gammaH_mean, r.gammaH_SD, 3),
                     "Hue method": r.hue_method})
    return pd.DataFrame(rows)


def Hmax_vs_gammaH(R):
    """Hmax (base colour) vs gammaH (colour change with angle), species means ± SD.

    Edge-method species have no comparable Hmax: their gammaH is shown in a side panel.
    """
    s = R["summary_by_species"].set_index("species")
    pk, ed = s[s.hue_method == "peak"], s[s.hue_method == "edge"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.4, 7), gridspec_kw=dict(width_ratios=[8, 1.3], wspace=.03)); tx = []
    for sp, d in pk.iterrows():
        a.errorbar(d.Hmax_nm_mean, d.gammaH_mean, xerr=d.Hmax_nm_SD, yerr=d.gammaH_SD, fmt="o", color=BLU, ms=8, capsize=2, lw=1.2, alpha=.9, zorder=3)
        tx.append(a.text(d.Hmax_nm_mean, d.gammaH_mean, name(sp), fontsize=10.5, style="italic"))
    adjust_text(tx, ax=a, arrowprops=dict(arrowstyle="-", color="grey", lw=.6))
    rho = spearmanr(pk.Hmax_nm_mean, pk.gammaH_mean)
    a.text(.03, .96, f"Spearman ρ = {rho.statistic:.2f}, p = {rho.pvalue:.2f}  (n = {len(pk)})", transform=a.transAxes, fontsize=12, color="#444", va="top")
    a.set_xlabel("$H_{max}$ — base colour, peak wavelength at span 0 (nm)"); a.set_ylabel("$\\gamma_H$ — hue change with angle (larger = more iridescent)")
    a.grid(alpha=.25); _clean(a)
    for i, (sp, d) in enumerate(ed.sort_index().iterrows()):
        x = (i - 1) * .45
        b.errorbar(x, d.gammaH_mean, yerr=d.gammaH_SD, fmt="s", mfc="white", color=RED, ms=7, capsize=2, lw=1.2)
        b.text(x, d.gammaH_mean + d.gammaH_SD + .004, name(sp), fontsize=9.5, style="italic", ha="center", va="bottom", rotation=90, color="#222")
    b.set_xlim(-.8, .8); b.set_ylim(a.get_ylim()); b.set_yticks([]); b.set_xticks([0]); b.set_xticklabels(["spectral edge\n($H_{max}$ n/a)"], fontsize=10); _clean(b); b.spines["left"].set_visible(False)
    return fig, dict(spearman_rho=rho.statistic, spearman_p=rho.pvalue, n=len(pk))
