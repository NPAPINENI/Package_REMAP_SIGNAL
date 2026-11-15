import argparse
import sys

import numpy as np
import pandas as pd

from .core import RemapSignal


def _regions(a):
    sep = "," if a.input.endswith(".csv") else "\t"
    df = pd.read_csv(a.input, sep=sep)
    with RemapSignal(a.species, a.bigwigs) as rs:
        res = rs.batch(df, stat=a.stat, one_based=a.one_based)
        if a.stat:
            res.to_csv(a.output, sep="\t", index=False)
        else:
            np.savez_compressed(a.output, tracks=np.array(rs.track_names),
                                **{f"r{i}": x for i, x in enumerate(res)})
    print(f"wrote {a.output}", file=sys.stderr)


def _build(a):
    from .pipeline import build_datasets
    build_datasets(a.species, a.bigwig, a.input_root, a.output_root,
                   input_filename=a.input_filename, tfs=a.tfs, seed=a.seed)


def _split(a):
    from .pipeline import make_splits
    make_splits(a.input, a.species, a.output_dir, bigwig=a.bigwig, name=a.name, seed=a.seed)


def main(argv=None):
    p = argparse.ArgumentParser(prog="remap-signal")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("regions", help="signal for regions in a TSV/CSV (chr,start,end[,sequence])")
    r.add_argument("--species", required=True, choices=["human", "mouse"])
    r.add_argument("--bigwigs", required=True, help="dir of .bw files, or one .bw file")
    r.add_argument("--input", required=True)
    r.add_argument("--output", required=True, help=".npz (per-base) or .tsv (summary)")
    r.add_argument("--stat", default=None, choices=["mean", "max", "sum", "min"])
    r.add_argument("--one-based", action="store_true")
    r.set_defaults(func=_regions)

    s = sub.add_parser("split", help="one file + species -> train/dev/test + sequence_*.tsv")
    s.add_argument("--input", required=True, help="TSV/CSV with chr,start,end,label,sequence")
    s.add_argument("--species", required=True, choices=["human", "mouse"])
    s.add_argument("--output-dir", required=True)
    s.add_argument("--bigwig", default=None, help="ReMap bigWig (default: $REMAP_HUMAN_BW / $REMAP_MOUSE_BW)")
    s.add_argument("--name", default=None)
    s.add_argument("--seed", type=int, default=42)
    s.set_defaults(func=_split)

    b = sub.add_parser("build", help="per-TF folders -> signal tables + train/dev/test splits")
    b.add_argument("--species", required=True, choices=["human", "mouse"])
    b.add_argument("--bigwig", required=True, help="single ReMap bigWig")
    b.add_argument("--input-root", required=True)
    b.add_argument("--output-root", required=True)
    b.add_argument("--input-filename", default="6mer_sequences.tsv")
    b.add_argument("--tfs", nargs="*", default=None, help="only these TF folders")
    b.add_argument("--seed", type=int, default=42)
    b.set_defaults(func=_build)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
