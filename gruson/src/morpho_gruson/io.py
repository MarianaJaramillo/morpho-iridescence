"""Read raw goniometer spreadsheets and parse their original column names.

Source headers may look like ``M. aega8 Brillo I12 C28.RAW8`` or
``aega8_matiz_I5C35.RAW8``. They contain a species, optional replicate number,
source sweep label, illumination angle I, collection angle C, and instrument
extension. The original source labels are deliberately retained in the parser;
parsed records use the English sweep types ``brightness`` and ``hue``.
White references are named ``Ref I<i> C<c>.RAW8`` / ``ref_I<i>C<c>.RAW8``.
"""
import re
import numpy as np
import openpyxl

SPECIES_ALIASES = {
    "rhetenro": "rhetenor", "rethenor": "rhetenor",
    "sulkowkyi": "sulkowskyi", "sulkowski": "sulkowskyi",
    "rhodopyeron": "rhodopteron",
    "godarti": "godartii", "godartti": "godartii",
    "helenoru": "helenor",
}
COLPAT = re.compile(
    r"^(?:M\.\s*)?([A-Za-zÑñ]+)\s*(\d*)(?:[_\s]*(brillo|matiz|brill|brilo))?_?\s*I?\s*(-?\d+)(?:\s*C\s*|\s+)(-?\d+)\s*\.?RAW8?$",
    re.IGNORECASE,
)
REFPAT = re.compile(r"^Ref\s*I\s*(-?\d+)\s*C\s*(-?\d+)\.RAW8$", re.IGNORECASE)
REFPAT2 = re.compile(r"^ref_I(-?\d+)C(-?\d+)\.RAW8$", re.IGNORECASE)


def load_sheet(path):
    """Return (wavelength[nm], {column name: intensity array}) from the first sheet."""
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    header = rows[0]
    start = 1
    for r in range(1, min(len(rows), 29)):
        if isinstance(rows[r][0], (int, float)):
            start = r
            break
    body = rows[start:]
    wl = np.array([r[0] for r in body], dtype=float)
    data = {}
    for j, h in enumerate(header):
        if j == 0 or not h:
            continue
        data[str(h).strip()] = np.array([r[j] if j < len(r) else None for r in body], dtype=float)
    return wl, data


def normalise_species(name):
    n = name.strip().lower()
    return SPECIES_ALIASES.get(n, n)


def parse_columns(names):
    """Split column names into sweep columns and white references.

    Returns (records, references) where
      records    : {(species, replicate): [(kind, I, C, column_name), ...]}  kind in {'brightness','hue'}
      references : [(I, C, column_name), ...]
    Unrecognised columns are ignored.
    """
    records, refs = {}, []
    for h in names:
        m = REFPAT.match(h) or REFPAT2.match(h)
        if m:
            refs.append((int(m.group(1)), int(m.group(2)), h))
            continue
        if re.search(r"blanc|white", h, re.I) or h.lower().strip() == "aana.raw8":
            continue
        m = COLPAT.match(h)
        if not m:
            continue
        sp_raw, rep, kind, I, C = m.groups()
        kind = "brightness" if kind and kind.lower().startswith("bril") else "hue"
        records.setdefault((normalise_species(sp_raw), rep), []).append((kind, int(I), int(C), h))
    return records, refs
