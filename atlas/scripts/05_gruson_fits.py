"""Gruson et al. (2019) model fits -> data/gruson/gruson_fits.json for the viewer.

Usage:
    python 05_gruson_fits.py <Gruson results folder> <data folder>

<Gruson results folder> holds the five tables, either in the legacy source folders parametros/ and ajustes/ or all in one folder (summary_by_species.csv, sweeps_brightness.csv,
sweeps_hue.csv, points_brightness.csv, points_hue.csv). A copy is in data/gruson/input/.

Model (Gruson et al. 2019, Interface Focus 9):
    brightness  B = Bmax * exp(-(x - t)^2 / (2 gammaB^2)),   x = bisector (I - C)/2, sweeps at constant span I + C
    hue         H = Hmax * cos(gammaH * span),                span = I + C (degrees), sweeps at constant I - C
Every fitted curve is recomputed from its parameters and checked against y_fitted.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

if len(sys.argv) != 3:
    sys.exit(__doc__)
D = os.path.join(sys.argv[1], '')
OUT = os.path.join(sys.argv[2], 'gruson')
os.makedirs(OUT, exist_ok=True)
def rd(f):
    # Keep these legacy source-folder names for compatibility; the pipeline outputs and documentation use English names.
    path = D + f if os.path.exists(D + f) else D + os.path.basename(f)
    return pd.read_csv(path, encoding='utf-8-sig')
S=rd('parametros/summary_by_species.csv');SB=rd('parametros/sweeps_brightness.csv');SH=rd('parametros/sweeps_hue.csv')
PB=rd('ajustes/points_brightness.csv');PH=rd('ajustes/points_hue.csv')
MAP={'godartii':'godarti'}   # accepted spelling used in the rest of the data
def sweep_key(df): return df.session.astype(str)+' #'+df.sweep.fillna('').astype(str).str.replace(r'\.0$','',regex=True)
for d in (SB,SH,PB,PH): d['key']=sweep_key(d).str.replace(' #$','',regex=True)
SB['c']=SB.constant_span_I_plus_C.round().astype(int); PB['c']=(PB.I_deg+PB.C_deg).astype(int)
SH['c']=SH.constant_I_minus_C.round().astype(int); PH['c']=(PH.I_deg-PH.C_deg).astype(int)
for d in (SB,SH,PB,PH): d['key2']=d.key+'|'+d.c.astype(str)
assert not SB.duplicated(['species','key2']).any() and not SH.duplicated(['species','key2']).any()
def assign(P,S):
    # each point goes to the sweep of the same session whose constant value is closest (some points are 1 degree off)
    out=[]
    for i,r in P.iterrows():
        c=S[(S.species==r.species)&(S.key==r.key)]
        out.append(c.key2.iloc[(c.c-r.c).abs().argmin()])
    P['key2']=out
assign(PB,SB); assign(PH,SH)
# chronological order of sessions
order=['31 ago-9 sep (S1)','2 sep','9 sep','25 sep','30 sep','7 oct','9 oct','15 oct','19 oct','21 oct']
def sess_rank(s): 
    for i,o in enumerate(order):
        if str(s).startswith(o): return i
    return 99
out={}
for sp in S.species:
    s=S[S.species==sp].iloc[0]
    def nz(x,d=3): return None if pd.isna(x) else round(float(x),d)
    sb=SB[SB.species==sp].copy(); sh=SH[SH.species==sp].copy()
    sb['r']=sb.session.map(sess_rank); sh['r']=sh.session.map(sess_rank)
    sb=sb.sort_values(['r','sweep']); sh=sh.sort_values(['r','sweep'])
    B=[];H=[]
    for _,r in sb.iterrows():
        p=PB[(PB.species==sp)&(PB.key2==r.key2)]
        B.append(dict(k=r.key,k2=r.key2,n=int(r.n_points),R2=round(r.R2,3),Bmax=round(r.Bmax,4),t=round(r.t_deg,3),g=round(r.gammaB_deg,3),span=nz(r.constant_span_I_plus_C,1),
            pts=[[int(a),int(b),round(float(x),2),round(float(y),4)] for a,b,x,y in zip(p.I_deg,p.C_deg,p.x_deg,p.y_observed)]))
    for _,r in sh.iterrows():
        p=PH[(PH.species==sp)&(PH.key2==r.key2)]
        H.append(dict(k=r.key,k2=r.key2,n=int(r.n_points),R2=round(r.R2,3),Hmax=round(r.Hmax_nm,2),g=round(r.gammaH,4),d=nz(r.constant_I_minus_C,1),note=None if pd.isna(r.note) else str(r.note),
            pts=[[int(a),int(b),round(float(x),2),round(float(y),2)] for a,b,x,y in zip(p.I_deg,p.C_deg,p.x_deg,p.y_observed)]))
    out[MAP.get(sp,sp)]=dict(name=sp,method=s.hue_method,nB=int(s.n_brightness),nH=int(s.n_hue),
        Bmax=[nz(s.Bmax_mean,4),nz(s.Bmax_SD,4)],t=[nz(s.t_deg_mean,3),nz(s.t_deg_SD,3)],gB=[nz(s.gammaB_deg_mean,3),nz(s.gammaB_deg_SD,3)],
        Hmax=[nz(s.Hmax_nm_mean,2),nz(s.Hmax_nm_SD,2)],gH=[nz(s.gammaH_mean,4),nz(s.gammaH_SD,4)],B=B,H=H)
    # check fitted curve reproduces y_fitted
    for r in B:
        p=PB[(PB.species==sp)&(PB.key2==r['k2'])]
        f=r['Bmax']*np.exp(-(p.x_deg-r['t'])**2/(2*r['g']**2)); assert np.allclose(f,p.y_fitted,atol=2e-3),(sp,r['k'],np.abs(f-p.y_fitted).max())
    for r in H:
        p=PH[(PH.species==sp)&(PH.key2==r['k2'])]
        f=r['Hmax']*np.cos(np.deg2rad(r['g']*p.x_deg))
        if not np.allclose(f,p.y_fitted,atol=0.05): print('hue model mismatch',sp,r['k'],float(np.abs(f-p.y_fitted).max()))
json.dump(out,open(os.path.join(OUT,'gruson_fits.json'),'w'),separators=(',',':'))
print(len(out),'species; sweeps B',sum(len(v['B']) for v in out.values()),'H',sum(len(v['H']) for v in out.values()))
print('max sweeps per species',max(max(len(v['B']),len(v['H'])) for v in out.values()))
