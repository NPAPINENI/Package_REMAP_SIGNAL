"""Per-TF dataset builder: attach ReMap signal to sequences and make train/dev/test splits.

Expected layout:  input_root/<TF>/<input_filename>   (TSV with chr,start,end,label,sequence)
Writes:           output_root/<TF>/Remap_signal_<TF>.tsv
                  output_root/<TF>/{train,dev,test}.tsv           (Sequence, Label)
                  output_root/<TF>/sequence_{train,dev,test}.tsv  (nucleotide_1..N signal)
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from .core import RemapSignal

REQUIRED = ["chr", "start", "end", "label", "sequence"]


def _fit(sig: np.ndarray, seq_len: int) -> np.ndarray:
    """Pad with zeros / truncate so the signal has exactly seq_len values."""
    out = np.zeros(seq_len, dtype=np.float32)
    n = min(len(sig), seq_len)
    out[:n] = sig[:n]
    return out


def signal_for_dataframe(rs: RemapSignal, df: pd.DataFrame, desc: str = "signal",
                         track: int = 0, progress: bool = True) -> pd.DataFrame:
    """Return DataFrame of nucleotide_1..N signal columns, one row per row of df.

    df needs chr, start, end, sequence. Signal length follows len(sequence).
    """
    sigs = []
    it = df.itertuples(index=False)
    if progress:
        it = tqdm(it, total=len(df), desc=desc)
    for r in it:
        L = len(r.sequence)
        try:
            s = rs.get_signal(r.chr, int(r.start), int(r.end), strict_length=False)[track]
            s = _fit(s, L)
        except Exception as e:  # keep going, like the original script
            print(f"Error {r.chr}:{r.start}-{r.end}: {e}")
            s = np.zeros(L, dtype=np.float32)
        sigs.append(s)

    width = max(len(s) for s in sigs) if sigs else 0
    mat = np.zeros((len(sigs), width), dtype=np.float32)
    for i, s in enumerate(sigs):
        mat[i, : len(s)] = s
    return pd.DataFrame(mat, columns=[f"nucleotide_{i+1}" for i in range(width)])


def split_70_15_15(final_df: pd.DataFrame, seed: int = 42):
    """Stratified 70/15/15 split on 'Label' (falls back to unstratified if classes too small)."""
    can = final_df["Label"].nunique() > 1 and final_df["Label"].value_counts().min() >= 2
    train, temp = train_test_split(
        final_df, test_size=0.30, stratify=final_df["Label"] if can else None,
        random_state=seed)
    can_t = can and temp["Label"].value_counts().min() >= 2
    dev, test = train_test_split(
        temp, test_size=0.50, stratify=temp["Label"] if can_t else None,
        random_state=seed)
    return train, dev, test


def process_file(input_file, rs: RemapSignal, out_dir, name: Optional[str] = None,
                 seed: int = 42, track: int = 0) -> dict:
    """One input table -> Remap_signal_<name>.tsv + train/dev/test + sequence_train/dev/test."""
    input_file = Path(input_file)
    name = name or input_file.stem
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {name} ===")

    sep = "," if input_file.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(input_file, sep=sep)
    df = df.rename(columns={"chrom": "chr", "seq": "sequence"})  # accept common aliases
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"[{name}] missing columns: {missing} (need {REQUIRED})")

    df = df.dropna(subset=REQUIRED).copy()
    df["sequence"] = df["sequence"].astype(str).str.upper()
    df["seq_len"] = df["sequence"].str.len()
    mism = (df["seq_len"] != (df["end"] - df["start"])).sum()
    if mism:
        print(f"[{name}] Warning: {mism} sequences differ in length from end-start")

    sig_df = signal_for_dataframe(rs, df, desc=f"{name}: signal", track=track)
    nuc_cols = list(sig_df.columns)
    seq_df = pd.DataFrame({"Sequence": df["sequence"].values, "Label": df["label"].values})
    final_df = pd.concat([seq_df, sig_df], axis=1)
    final_df.to_csv(out_dir / f"Remap_signal_{name}.tsv", sep="\t", index=False)

    train, dev, test = split_70_15_15(final_df, seed)
    for split, d in [("train", train), ("dev", dev), ("test", test)]:
        d[["Sequence", "Label"]].to_csv(out_dir / f"{split}.tsv", sep="\t", index=False)
        d[nuc_cols].to_csv(out_dir / f"sequence_{split}.tsv", sep="\t", index=False)

    print(f"[{name}] total {len(final_df):,} | train {len(train):,} | "
          f"dev {len(dev):,} | test {len(test):,}")
    return {"all": final_df, "train": train, "dev": dev, "test": test}


def process_tf(tf_dir: Path, rs: RemapSignal, output_root: Path,
               input_filename: str = "6mer_sequences.tsv", seed: int = 42,
               track: int = 0):
    tf_dir = Path(tf_dir)
    return process_file(tf_dir / input_filename, rs, Path(output_root) / tf_dir.name,
                        name=tf_dir.name, seed=seed, track=track)["all"]


def make_splits(input_file, species: str, output_dir, bigwig: Optional[str] = None,
                name: Optional[str] = None, seed: int = 42) -> dict:
    """Main entry point: any file (chr,start,end,label,sequence) + species -> splits + ReMap signal.

    bigwig=None falls back to the REMAP_HUMAN_BW / REMAP_MOUSE_BW environment variable. Returns dict of DataFrames: all/train/dev/test.
    """
    from .multi import resolve_bigwig
    from .core import _norm_species
    sp = _norm_species(species)
    with RemapSignal(sp, resolve_bigwig(sp, bigwig)) as rs:
        return process_file(input_file, rs, output_dir, name=name, seed=seed)


def build_datasets(
    species: str,
    bigwig: str,
    input_root: Union[str, Path],
    output_root: Union[str, Path],
    input_filename: str = "6mer_sequences.tsv",
    tfs: Optional[List[str]] = None,
    seed: int = 42,
) -> None:
    """Process every input_root/<TF>/<input_filename> for the given species + ReMap bigWig."""
    input_root, output_root = Path(input_root), Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    tf_dirs = sorted(p.parent for p in input_root.glob(f"*/{input_filename}"))
    if tfs:
        tf_dirs = [d for d in tf_dirs if d.name in set(tfs)]
    if not tf_dirs:
        raise FileNotFoundError(f"No '{input_filename}' found under {input_root}/<TF>/")
    print(f"Found {len(tf_dirs)} TF folders: {[d.name for d in tf_dirs]}")

    with RemapSignal(species, bigwig) as rs:
        for d in tf_dirs:
            process_tf(d, rs, output_root, input_filename, seed)
    print("\nAll TF folders processed.")
