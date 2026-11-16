"""Find (or build) the ReMap bigWig for a species and query human + mouse together."""
from __future__ import annotations

import os
from typing import Optional

from .core import RemapSignal, _norm_species

LEGACY_ENV = {"hg38": "REMAP_HUMAN_BW", "mm39": "REMAP_MOUSE_BW"}


def resolve_bigwig(species: str, bigwig: Optional[str] = None) -> str:
    """Explicit path > env var (REMAP_<ASSEMBLY>_BW, e.g. REMAP_HUMAN_BW) > cached > auto-build."""
    from .download import cached_path, download
    asm = _norm_species(species)
    path = bigwig
    if not path:
        for name in (f"REMAP_{asm.upper()}_BW", LEGACY_ENV.get(asm)):
            if name and os.environ.get(name):
                path = os.environ[name]
                break
    if path:
        if not os.path.exists(path):
            raise FileNotFoundError(f"ReMap bigWig not found: {path}")
        return path
    cached = cached_path(asm)
    return cached if os.path.exists(cached) else download(asm)


class Remap:
    """Human (hg38) and mouse (mm39) together; pick the species per call.

    Files are found/built lazily on first use of each species.
    """

    def __init__(self, human: Optional[str] = None, mouse: Optional[str] = None, **others):
        self._given = {"hg38": human, "mm39": mouse, **{_norm_species(k): v for k, v in others.items()}}
        self._rs = {}

    def __getitem__(self, species: str) -> RemapSignal:
        asm = _norm_species(species)
        if asm not in self._rs:
            self._rs[asm] = RemapSignal(asm, resolve_bigwig(asm, self._given.get(asm)))
        return self._rs[asm]

    human = property(lambda self: self["human"])
    mouse = property(lambda self: self["mouse"])

    def get_signal(self, species, chrom, start, end, sequence=None, **kw):
        return self[species].get_signal(chrom, start, end, sequence, **kw)

    def get_signal_df(self, species, chrom, start, end, sequence=None, **kw):
        return self[species].get_signal_df(chrom, start, end, sequence, **kw)

    def summarize(self, species, chrom, start, end, sequence=None, **kw):
        return self[species].summarize(chrom, start, end, sequence, **kw)

    def close(self):
        for r in self._rs.values():
            r.close()
        self._rs.clear()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def load(human: Optional[str] = None, mouse: Optional[str] = None, **others) -> Remap:
    """Shortcut: remap_signal.load() -> Remap ready for human and mouse."""
    return Remap(human, mouse, **others)
