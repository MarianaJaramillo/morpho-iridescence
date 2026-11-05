"""Collects every correction and every detected problem, with the file and column it came from."""
import pandas as pd

_entries = []


def log(species_id, variant, file, kind, detail, column='-'):
    _entries.append(dict(species_id=species_id, variant=variant, file=file, column=column, kind=kind, detail=detail))


def to_frame():
    return pd.DataFrame(_entries, columns=['species_id', 'variant', 'file', 'column', 'kind', 'detail'])
