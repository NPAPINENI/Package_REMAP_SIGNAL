"""Fetch ReMap bigWig signal for (chr, start, end, sequence) in human or mouse.

Coordinates are 0-based, half-open (BED / bigWig convention) unless
``one_based=True`` is passed (then start/end are 1-based inclusive).
"""
from __future__ import annotations

import glob
import os
from typing import Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
import pyBigWig

ASSEMBLIES = ("hg38", "hg19", "mm10", "mm39", "dm6")
SPECIES_ALIASES = {
    "human": "hg38", "hs": "hg38", "homo_sapiens": "hg38",
    "mouse": "mm39", "mm": "mm39", "mus_musculus": "mm39",
    "fly": "dm6", "drosophila": "dm6",
    **{a: a for a in ASSEMBLIES},
}


def _norm_species(species: str) -> str:
    """'human' -> 'hg38', 'mouse' -> 'mm39'; assembly names (hg19, mm10, dm6, ...) pass through."""
    try:
        return SPECIES_ALIASES[species.lower()]
    except KeyError:
        raise ValueError(
            f"Unknown species/assembly {species!r}. Use human, mouse, fly or one of {ASSEMBLIES}")


class RemapSignal:
    """Query ReMap bigWig tracks for one species.

    Parameters
    ----------
    species : "human" | "mouse" (aliases: hg38, mm10, mm39, ...)
    bigwigs : one of
        - directory containing ``*.bw`` / ``*.bigWig`` files (track name = file stem)
        - dict {track_name: path}
        - list of paths
        - single path to one merged bigWig
    fill_nan : value used where the bigWig has no data (default 0.0)
    """

    def __init__(
        self,
        species: str,
        bigwigs: Union[str, Dict[str, str], Sequence[str]],
        fill_nan: float = 0.0,
    ):
        try:
            self.species = _norm_species(species)
        except ValueError:          # any other species is fine when you supply the bigWig
            self.species = species
        self.fill_nan = fill_nan
        self.paths = self._resolve(bigwigs)
        if not self.paths:
            raise FileNotFoundError("No bigWig files found")
        self._handles: Dict[str, pyBigWig.pyBigWig] = {}
        for name, p in self.paths.items():
            bw = pyBigWig.open(p)
            if bw is None:
                raise OSError(f"Cannot open bigWig: {p}")
            self._handles[name] = bw
        self.track_names: List[str] = list(self._handles)

    # ------------------------------------------------------------------ setup
    @staticmethod
    def _resolve(bigwigs) -> Dict[str, str]:
        if isinstance(bigwigs, dict):
            return dict(bigwigs)
        if isinstance(bigwigs, str):
            if os.path.isdir(bigwigs):
                files = []
                for ext in ("*.bw", "*.bigWig", "*.bigwig"):
                    files += glob.glob(os.path.join(bigwigs, ext))
                bigwigs = sorted(files)
            else:
                bigwigs = [bigwigs]
        return {os.path.splitext(os.path.basename(p))[0]: p for p in bigwigs}

    def close(self):
        for bw in self._handles.values():
            bw.close()
        self._handles.clear()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    # ------------------------------------------------------------- chrom names
    def _chrom_for(self, bw, chrom: str) -> str:
        chroms = bw.chroms()
        if chrom in chroms:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in chroms:
            return alt
        if chrom in ("M", "chrM", "MT", "chrMT") :
            for c in ("chrM", "M", "MT", "chrMT"):
                if c in chroms:
                    return c
        raise KeyError(f"{chrom!r} not in bigWig {sorted(chroms)[:5]}...")

    # ------------------------------------------------------------------- query
    def get_signal(
        self,
        chrom: str,
        start: int,
        end: int,
        sequence: Optional[str] = None,
        one_based: bool = False,
        strict_length: bool = True,
    ) -> np.ndarray:
        """Return array of shape (n_tracks, end-start): per-base ReMap signal."""
        if one_based:
            start -= 1
        length = end - start
        if length <= 0:
            raise ValueError(f"end must be > start ({chrom}:{start}-{end})")
        if sequence is not None and len(sequence) != length:
            msg = (f"sequence length {len(sequence)} != end-start {length} "
                   f"for {chrom}:{start}-{end}")
            if strict_length:
                raise ValueError(msg)
            import warnings
            warnings.warn(msg)

        out = np.full((len(self._handles), length), self.fill_nan, dtype=np.float32)
        for i, bw in enumerate(self._handles.values()):
            try:
                c = self._chrom_for(bw, chrom)
            except KeyError:
                continue  # chromosome absent -> all fill value
            size = bw.chroms(c)
            s, e = max(start, 0), min(end, size)
            if e <= s:
                continue
            vals = np.asarray(bw.values(c, s, e, numpy=True), dtype=np.float32)
            vals[np.isnan(vals)] = self.fill_nan
            out[i, s - start : e - start] = vals
        return out

    def get_signal_df(self, chrom, start, end, sequence=None, **kw) -> pd.DataFrame:
        """Per-base table: chrom, pos, base (if sequence given), one column per track."""
        arr = self.get_signal(chrom, start, end, sequence, **kw)
        s = start - 1 if kw.get("one_based") else start
        df = pd.DataFrame(arr.T, columns=self.track_names)
        df.insert(0, "pos", np.arange(s, s + arr.shape[1]))
        df.insert(0, "chrom", chrom)
        if sequence is not None:
            df.insert(2, "base", list(sequence))
        return df

    def summarize(self, chrom, start, end, sequence=None, stat="mean", **kw) -> pd.Series:
        """One value per track over the region (mean/max/sum/min)."""
        arr = self.get_signal(chrom, start, end, sequence, **kw)
        f = {"mean": np.mean, "max": np.max, "sum": np.sum, "min": np.min}[stat]
        return pd.Series(f(arr, axis=1), index=self.track_names)

    def batch(
        self,
        regions: Union[pd.DataFrame, Iterable[tuple]],
        stat: Optional[str] = None,
        **kw,
    ):
        """Run many regions.

        regions: DataFrame with columns chr/chrom, start, end, [sequence]
                 or iterable of (chrom, start, end[, sequence]).
        stat=None  -> list of (n_tracks, L) arrays
        stat="mean"/"max"/"sum"/"min" -> DataFrame (rows=regions, cols=tracks)
        """
        if isinstance(regions, pd.DataFrame):
            df = regions.rename(columns={"chr": "chrom", "seq": "sequence"})
            if "sequence" not in df:
                df["sequence"] = None
            rows = list(zip(df["chrom"], df["start"], df["end"], df["sequence"]))
        else:
            rows = [tuple(r) + (None,) * (4 - len(r)) for r in regions]
        if stat is None:
            return [self.get_signal(c, s, e, q, **kw) for c, s, e, q in rows]
        res = [self.summarize(c, s, e, q, stat=stat, **kw) for c, s, e, q in rows]
        out = pd.DataFrame(res)
        out.insert(0, "end", [r[2] for r in rows])
        out.insert(0, "start", [r[1] for r in rows])
        out.insert(0, "chrom", [r[0] for r in rows])
        return out.reset_index(drop=True)
