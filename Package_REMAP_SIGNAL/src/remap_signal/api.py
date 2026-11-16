"""Simple one-line helpers: remap_signal.get_signal(...) and remap_signal.add_signal(...)."""
from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd

from .core import RemapSignal
from .multi import resolve_bigwig

_open = {}


def _reader(species: str, bigwig: Optional[str]) -> RemapSignal:
    key = (species.lower(), bigwig)
    if key not in _open:
        if bigwig:
            if not os.path.exists(bigwig):
                raise FileNotFoundError(f"bigWig not found: {bigwig}")
            path = bigwig
        else:
            path = resolve_bigwig(species)
        _open[key] = RemapSignal(species, path)
    return _open[key]


def get_signal(species: str, chrom: str, start: int, end: int,
               bigwig: Optional[str] = None) -> np.ndarray:
    """ReMap signal for one region: a 1-D array with one value per base (end - start).

    species: "human", "mouse", or any name if you give bigwig="/path/file.bw".
    Coordinates are 0-based, end exclusive.
    """
    arr = _reader(species, bigwig).get_signal(chrom, int(start), int(end))
    return arr[0] if arr.shape[0] == 1 else arr


def add_signal(df: pd.DataFrame, species: str, bigwig: Optional[str] = None) -> pd.DataFrame:
    """Return df with the signal appended as columns nucleotide_1 ... nucleotide_N.

    df needs the columns chr, start, end. Shorter regions are padded with 0.
    """
    d = df.rename(columns={"chrom": "chr"})
    r = _reader(species, bigwig)
    sigs = [r.get_signal(c, int(s), int(e))[0] for c, s, e in zip(d["chr"], d["start"], d["end"])]
    width = max((len(x) for x in sigs), default=0)
    mat = np.zeros((len(sigs), width), dtype=np.float32)
    for i, x in enumerate(sigs):
        mat[i, : len(x)] = x
    cols = [f"nucleotide_{i + 1}" for i in range(width)]
    sig = pd.DataFrame(mat, columns=cols, index=df.index)
    return pd.concat([df, sig], axis=1)
