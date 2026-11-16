"""Get the ReMap 2022 density bigWig for a species (downloaded once, then cached).

Default ("all"): the UCSC ReMap density track
  https://hgdownload.soe.ucsc.edu/gbdb/<assembly>/reMap/reMapDensity2022.bw
for hg38 (human) and mm39 (mouse). It is downloaded directly and checked to be a valid bigWig.
Other assemblies (mm10, hg19, dm6) are built from the ReMap website BED file instead:
  kind "all" = summit +/- 25 bp, "nr" / "crm" = full BED interval, counted per base.
Choose with download(kind=...) or the REMAP_KIND environment variable.
Cache folder: $REMAP_CACHE_DIR, default ~/.cache/remap_signal.
"""
from __future__ import annotations

import os
import shutil
import urllib.request
from typing import Optional

import numpy as np
import pandas as pd
import pyBigWig
from tqdm import tqdm

from .core import _norm_species

# ready-made UCSC bigWigs (md5 unknown -> validity is checked after download)
UCSC_URL = "https://hgdownload.soe.ucsc.edu/gbdb/{asm}/reMap/reMapDensity2022.bw"
DIRECT = {
    "hg38": {"url": UCSC_URL.format(asm="hg38"), "md5": None},
    "mm39": {"url": UCSC_URL.format(asm="mm39"), "md5": None},
}
FLANK = 25   # summit +/- 25 bp  (the "w50" recipe)

BED_URL = ("https://remap.univ-amu.fr/storage/remap2022/{asm}/MACS2/"
           "remap2022_{kind}_macs2_{asm}_v1_0.bed.gz")
KINDS = ("crm", "nr", "all")


def default_kind() -> str:
    return os.environ.get("REMAP_KIND", "all")


def cache_dir() -> str:
    d = os.environ.get("REMAP_CACHE_DIR") or os.path.join(
        os.path.expanduser("~"), ".cache", "remap_signal")
    os.makedirs(d, exist_ok=True)
    return d


def cached_path(species: str, kind: Optional[str] = None) -> str:
    asm, kind = _norm_species(species), kind or default_kind()
    return os.path.join(cache_dir(), f"{asm}_remap2022_{kind}_density.bw")


def _fetch(url: str, dest: str, desc: str):
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(url) as r, open(tmp, "wb") as out:
            total = int(r.headers.get("Content-Length") or 0) or None
            with tqdm(total=total, unit="B", unit_scale=True, desc=desc) as bar:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    bar.update(len(chunk))
        shutil.move(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _md5(path: str) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bed_to_density_bigwig(bed_gz: str, out_bw: str, summit: bool = False, flank: int = FLANK):
    """Write a bigWig: value = number of regions covering each base.

    summit=False: regions are BED columns 2-3.
    summit=True : regions are the peak summit (BED columns 7-8) extended by `flank` on each side.
    """
    starts, ends = {}, {}
    reader = pd.read_csv(bed_gz, sep="\t", header=None, usecols=[0, 6, 7] if summit else [0, 1, 2],
                         names=["c", "s", "e"], dtype={"c": str, "s": np.int64, "e": np.int64},
                         comment="#", chunksize=5_000_000)
    for ch in reader:
        if summit:
            ch["s"] = (ch["s"] - flank).clip(lower=0)
            ch["e"] = ch["e"] + flank
        ch = ch[ch["c"] != "chrEBV"]
        for c, g in ch.groupby("c", sort=False):
            starts.setdefault(c, []).append(g["s"].values)
            ends.setdefault(c, []).append(g["e"].values)

    chroms = sorted(starts)
    tmp = out_bw + ".part"
    bw = pyBigWig.open(tmp, "w")
    sizes = {c: int(max(a.max() for a in ends[c])) + 1 for c in chroms}
    bw.addHeader([(c, sizes[c]) for c in chroms])
    for c in chroms:
        s, e = np.concatenate(starts.pop(c)), np.concatenate(ends.pop(c))
        pos = np.concatenate([s, e])
        delta = np.concatenate([np.ones(len(s)), -np.ones(len(e))])
        u, inv = np.unique(pos, return_inverse=True)
        cov = np.cumsum(np.bincount(inv, weights=delta))
        keep = cov[:-1] > 0
        if keep.any():
            st = u[:-1][keep].tolist()
            bw.addEntries([c] * len(st), st, ends=u[1:][keep].tolist(),
                          values=cov[:-1][keep].astype(float).tolist())
    bw.close()
    shutil.move(tmp, out_bw)


def download(species: str, kind: Optional[str] = None, force: bool = False,
             url: Optional[str] = None) -> str:
    """Download ReMap peaks (once) and build the density bigWig. Returns its local path."""
    asm, kind = _norm_species(species), kind or default_kind()
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    dest = cached_path(asm, kind)
    if os.path.exists(dest) and not force:
        return dest

    if url is None and kind == "all" and asm in DIRECT:      # ready-made UCSC bigWig
        d = DIRECT[asm]
        print(f"Downloading ReMap 2022 {asm} density bigWig from {d['url']} (one time, large file)")
        new = dest + ".dl"
        _fetch(d["url"], new, asm)
        bad = d["md5"] and _md5(new) != d["md5"]
        if not bad:
            try:
                h = pyBigWig.open(new)
                bad = h is None or not h.isBigWig()
                if h is not None:
                    h.close()
            except Exception:
                bad = True
        if bad:
            os.remove(new)                      # keep any earlier good cached file
            raise OSError(f"Downloaded {asm} file is not a valid bigWig (check the URL): {d['url']}")
        shutil.move(new, dest)
        return dest

    url = url or BED_URL.format(asm=asm, kind=kind)
    bed = dest.replace("_density.bw", ".bed.gz")
    print(f"Downloading ReMap {asm} ({kind}) peaks from {url}")
    if kind == "all":
        print("  note: the 'all' file is large (several GB) and building needs a lot of memory")
    _fetch(url, bed, asm)
    print("Building bigWig (one time)...")
    bed_to_density_bigwig(bed, dest, summit=(kind == "all"))
    os.remove(bed)
    print(f"  -> {dest}")
    return dest
