"""Load human and mouse ReMap bigWigs together and query by species."""
from __future__ import annotations

import os
from typing import Optional

from .core import RemapSignal, _norm_species

ENV_VARS = {"human": "REMAP_HUMAN_BW", "mouse": "REMAP_MOUSE_BW"}


def resolve_bigwig(species: str, bigwig: Optional[str] = None) -> str:
    """Explicit path > environment variable (REMAP_HUMAN_BW / REMAP_MOUSE_BW)."""
    sp = _norm_species(species)
    path = bigwig or os.environ.get(ENV_VARS[sp])
    if not path:
        raise ValueError(
            f"No ReMap bigWig for {sp}. Pass bigwig='/path/file.bw' or set "
            f"the environment variable {ENV_VARS[sp]}."
        )
    if not os.path.exists(path):
        raise FileNotFoundError(f"ReMap bigWig not found: {path}")
    return path


class Remap:
    """Holds a human and a mouse RemapSignal; pick the species per call.

    Paths: argument > env var (REMAP_HUMAN_BW / REMAP_MOUSE_BW).
    Species are opened lazily, so a missing mouse file doesn't break human-only use.
    """

    def __init__(self, human: Optional[str] = None, mouse: Optional[str] = None):
        self._given = {"human": human, "mouse": mouse}
        self._rs = {}

    def __getitem__(self, species: str) -> RemapSignal:
        sp = _norm_species(species)
        if sp not in self._rs:
            self._rs[sp] = RemapSignal(sp, resolve_bigwig(sp, self._given[sp]))
        return self._rs[sp]

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


def load(human: Optional[str] = None, mouse: Optional[str] = None) -> Remap:
    """Shortcut: remap_signal.load() -> Remap with both species ready."""
    return Remap(human, mouse)
