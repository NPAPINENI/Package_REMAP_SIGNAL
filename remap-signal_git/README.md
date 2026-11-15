# remap-signal

Give it a table of genomic regions with sequences and a species (human or mouse).
It returns the per-base **ReMap** signal for every sequence and writes stratified
**train / dev / test** splits.

## Install

```bash
conda install -c bioconda -c conda-forge pybigwig   # easiest way to get pyBigWig
pip install git+https://github.com/YOUR_USERNAME/remap-signal.git
```

Python >= 3.8.

## You need the ReMap bigWigs

Download the ReMap 2022 density bigWig for your genome (human hg38, mouse mm10) from the
ReMap website. Then either pass the path each time or set it once:

```bash
export REMAP_HUMAN_BW=/path/to/human/reMapDensity2022.bw
export REMAP_MOUSE_BW=/path/to/mouse/reMapDensity2022.bw
```

Your coordinates must use the same assembly as the bigWig (0-based, half-open).

## Input file

TSV or CSV with columns `chr, start, end, label, sequence` (`chrom` and `seq` are accepted
aliases). Extra columns are ignored.

## Quick start

```python
import remap_signal as rs

rs.make_splits("my_file.tsv", species="human", output_dir="out_human/")
rs.make_splits("my_mouse.tsv", species="mouse", output_dir="out_mouse/",
               bigwig="/path/to/mouse.bw")      # optional override
```

or from the shell:

```bash
remap-signal split --input my_file.tsv --species human --output-dir out_human/
```

### Output (in `output_dir`)

| File | Content |
|---|---|
| `train.tsv`, `dev.tsv`, `test.tsv` | `Sequence`, `Label` |
| `sequence_train.tsv`, `sequence_dev.tsv`, `sequence_test.tsv` | ReMap signal, columns `nucleotide_1..N`, row-aligned with the files above |
| `Remap_signal_<name>.tsv` | everything before splitting |

Split: stratified 70/15/15 on `Label`, seed 42 (`seed=` to change).

## Other uses

```python
# signal for one region (rows = tracks, columns = bases)
remap = rs.load()                       # uses the env vars above
arr = remap.get_signal("human", "chr1", 1_000_000, 1_000_300, seq)
df  = remap.get_signal_df("mouse", "chr2", 5_000_000, 5_000_300, seq)
remap.summarize("human", "chr1", 1_000_000, 1_000_300, stat="max")

# many <TF>/6mer_sequences.tsv folders at once
rs.build_datasets("human", "/path/human.bw", input_root="data/", output_root="out/")
```

Chromosome names `1` and `chr1` are matched either way. Positions with no signal, or
outside the chromosome, are 0.

## Develop

```bash
git clone https://github.com/YOUR_USERNAME/remap-signal.git && cd remap-signal
pip install -e ".[dev]"
pytest
```

## License

MIT
