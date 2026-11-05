"""Step 1 - header repair.

Each raw column is named after its geometry, e.g. 'M. absoloni -10I 20C.RFL8'. The names were typed by hand and
contain typos. This step turns every name into (I, C) in degrees and logs each correction:

  * '+' inside a number            '-1+0I'  -> '-10I'
  * '+C'                            '20+C'   -> '20C'
  * letter O instead of zero        '1OC'    -> '10C'
  * missing space between I and C  '-10I20C' -> '-10I 20C'
  * missing letter C               '-10I 20' -> C = 20
  * I and C swapped in a whole file (the 5-degree sweep must be C)
  * C out of range (|C| > 90)       '-125C'  -> the last two digits, C = -25
  * C off the 5-degree grid ('4C', '-79C'): assigned to the empty slot of that incidence row whose text is most similar
"""
import difflib
import re

import numpy as np

from .changelog import log
from .config import C_GRID


def parse_header(name):
    """Return (I, C, fixes) from one raw spectrum name. I or C is None when it cannot be read."""
    s = re.sub(r'\.RFL8$', '', str(name).strip(), flags=re.I)
    fixes = []
    m = re.match(r'^(.*?)\s+(-?\d.*)$', s)
    r = m.group(2) if m else s
    if re.search(r'\d\+\d', r):
        r = re.sub(r'(\d)\+(\d)', r'\1\2', r); fixes.append('"+" sign inside the number')
    if '+C' in r:
        r = r.replace('+C', 'C'); fixes.append('"+C" read as "C"')
    if re.search(r'\dO[CI]\b', r):
        r = re.sub(r'(\d)O([CI])', r'\g<1>0\2', r); fixes.append('letter O instead of zero')
    if re.search(r'\d[IC]-?\d', r):
        r = re.sub(r'(\d[IC])(-?\d)', r'\1 \2', r); fixes.append('missing space between I and C')
    tokens = re.findall(r'(-?\d+)\s*([IC])?', r)
    d, bare = {}, []
    for v, k in tokens:
        if k:
            d[k] = float(v)
        else:
            bare.append(float(v))
    if 'I' in d and 'C' not in d and len(bare) == 1:
        d['C'] = bare[0]; fixes.append('missing letter C')
    return d.get('I'), d.get('C'), fixes


def repair_file_headers(names, species_id, variant, file_name):
    """Parse all names of one file. Returns a list of (I, C) and logs every fix."""
    parsed = [parse_header(n) for n in names]
    Is = {p[0] for p in parsed if p[0] is not None}
    Cs = {p[1] for p in parsed if p[1] is not None}
    swapped = len(Is) > len(Cs) and any(v % 10 for v in Is)
    if swapped:
        log(species_id, variant, file_name, 'I_C_swapped',
            'Names written with I and C the wrong way round (the 5-degree sweep was labelled I). Swapped back.', 'all')
    out = []
    for j, (n, (I, C, fixes)) in enumerate(zip(names, parsed)):
        if swapped:
            I, C = C, I
        for fx in fixes:
            log(species_id, variant, file_name, 'header_fixed', f'"{n}": {fx}', j + 1)
        if C is not None and abs(C) > 90:
            newC = float(str(int(C))[-2:]) * np.sign(C)
            log(species_id, variant, file_name, 'header_fixed', f'"{n}": C={C:g} out of range, read as {newC:g}', j + 1)
            C = newC
        out.append((I, C))
    return out


def assign_off_grid(meta, species_id, variant):
    """Move angles that are not on the 5-degree grid to the free slot of their row with the most similar text."""
    off = meta[(meta.C % 5 != 0)]
    for i, r in off.iterrows():
        row = meta[(meta.I == r.I)]
        holes = [c for c in C_GRID if c not in set(row.C) and abs(c - r.I) > 5]
        if holes:
            best = max(holes, key=lambda c: difflib.SequenceMatcher(None, f'{int(r.C)}', f'{int(c)}').ratio())
            log(species_id, variant, r.source_file, 'header_fixed',
                f'"{r.original_header}": C={r.C:g} is not on the grid; assigned to the empty slot C={best:g} of row I={r.I:g}',
                r.source_column)
            meta.loc[i, 'C'] = float(best)
    return meta
