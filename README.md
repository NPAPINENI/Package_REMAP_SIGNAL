# Package: REMAP_SIGNAL

### instructions to get this package to your environment
#### 1. Activate your environment
#### 2. conda install -c bioconda -c conda-forge pybigwig scikit-learn pandas numpy tqdm -y
#### 3. python -m pip install git+https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git


Motivation: To attach **ReMap** ChIP-seq signal to genomic sequences (human or mouse) and build
**train / dev / test** splits, in one call.

Result: Give it a table with `chr, start, end, label, sequence` and a species. You get back the
per-base ReMap signal for every sequence, aligned to the sequence, plus stratified
70/15/15 splits.

---

## 1. Requirements

- Python 3.8 or newer
- The ReMap bigWig for each species you use: human (hg38) and/or mouse (mm10).
  Download the ReMap 2022 density bigWigs from the ReMap website.
- Python packages (installed automatically except pyBigWig, see below):
  numpy, pandas, pyBigWig, scikit-learn, tqdm

## 2. Install

### 2.1 Create or activate an environment

```bash
conda create -n remap python=3.10 -y
conda activate remap
```

### 2.2 Install the dependencies with conda

pyBigWig is the one that most often fails with pip, so install it with conda:

```bash
conda install -c bioconda -c conda-forge pybigwig scikit-learn pandas numpy tqdm -y
```

### 2.3 Install the package from GitHub

```bash
python -m pip install git+https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git
```

Use `python -m pip`, not plain `pip`, so the package goes into the active environment.
On shared servers, plain `pip` is often the system one and fails with
`externally-managed-environment`.

If the package files sit in a subfolder of the repository, add the folder name:

```bash
python -m pip install "git+https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git#subdirectory=<folder>"
```

### 2.4 Check the install

```bash
python -c "import remap_signal; print(remap_signal.__file__)"
remap-signal --help
```

### 2.5 Update later

```bash
python -m pip install --upgrade --force-reinstall --no-deps git+https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git
```

## 3. Tell the package where the ReMap bigWigs are

Set two environment variables once (put them in `~/.bashrc`, then run `source ~/.bashrc`):

```bash
export REMAP_HUMAN_BW=/path/to/human/reMapDensity2022.bw
export REMAP_MOUSE_BW=/path/to/mouse/reMapDensity2022.bw
```

Or skip this and pass the path every time (`bigwig="..."` in Python, `--bigwig ...` on the
command line). An explicit path always overrides the environment variable.

Your coordinates must use the **same genome assembly** as the bigWig (human hg38, mouse
mm10). Coordinates are **0-based, half-open** (BED style).

## 4. Input file

A TSV or CSV with these columns:

| column | meaning |
|---|---|
| `chr` | chromosome, e.g. `chr1` (`1` also works) |
| `start` | start position (0-based) |
| `end` | end position (exclusive) |
| `label` | class label, used for the stratified split |
| `sequence` | DNA sequence, same length as `end - start` |

`chrom` and `seq` are accepted as alternative names for `chr` and `sequence`.
Any extra columns are ignored. Rows with missing values in these columns are dropped.

Example (tab-separated):

```
chr	start	end	label	sequence
chr1	1000000	1000300	1	ACGT...
chr2	5000000	5000300	0	TTGA...
```

## 5. Quick start

### Python

```python
import remap_signal as rs

# human
rs.make_splits("my_file.tsv", species="human", output_dir="out_human/")

# mouse
rs.make_splits("my_mouse_file.tsv", species="mouse", output_dir="out_mouse/")

# explicit bigWig, custom name and seed
rs.make_splits("my_file.tsv", species="human", output_dir="out/",
               bigwig="/path/to/human.bw", name="my_run", seed=42)
```

`make_splits` also returns a dictionary of DataFrames with the keys
`all`, `train`, `dev` and `test`.

### Command line

```bash
remap-signal split --input my_file.tsv --species human --output-dir out_human/
remap-signal split --input my_mouse.tsv --species mouse --output-dir out_mouse/ \
                   --bigwig /path/to/mouse.bw
```

Options: `--bigwig`, `--name`, `--seed`.

## 6. Output

Written to `output_dir`:

| file | content |
|---|---|
| `train.tsv`, `dev.tsv`, `test.tsv` | columns `Sequence`, `Label` |
| `sequence_train.tsv`, `sequence_dev.tsv`, `sequence_test.tsv` | ReMap signal, columns `nucleotide_1 ... nucleotide_N` |
| `Remap_signal_<name>.tsv` | everything before splitting |

Row *i* of `sequence_train.tsv` is the signal for row *i* of `train.tsv`
(same for dev and test).

The split is a stratified 70/15/15 on `Label` with seed 42. If a class is too small to
stratify, it falls back to a random split. Positions with no signal, or outside the
chromosome, are 0.

## 7. Other things you can do

### Signal for one region

```python
import remap_signal as rs

remap = rs.load()   # uses REMAP_HUMAN_BW / REMAP_MOUSE_BW, opens each species when first used

arr = remap.get_signal("human", "chr1", 1_000_000, 1_000_300, seq)
# numpy array, shape (n_tracks, 300): one row per bigWig, one column per base

df = remap.get_signal_df("mouse", "chr2", 5_000_000, 5_000_300, seq)
# per-base table: chrom, pos, base, signal

remap.summarize("human", "chr1", 1_000_000, 1_000_300, stat="max")
# one value per track: stat can be mean, max, sum or min
```

`seq` is optional. When given, its length must equal `end - start`.

### One species, your own bigWig(s)

```python
from remap_signal import RemapSignal

with RemapSignal("human", "/path/to/file.bw") as r:
    arr = r.get_signal("chr1", 1_000_000, 1_000_300)

# several tracks: a folder of .bw files, a list, or {name: path}
r = RemapSignal("mouse", {"Ctcf": "/p/Ctcf.bw", "Sox2": "/p/Sox2.bw"})
```

### Many regions at once

```python
res = r.batch(regions_df, stat="mean")   # DataFrame: one row per region, one column per track
arrs = r.batch(regions_df)               # list of per-base arrays
```

`regions_df` needs the columns `chr, start, end` and optionally `sequence`.

### Command line: signal for regions

```bash
remap-signal regions --species human --bigwigs /path/to/file.bw \
                     --input regions.tsv --output out.tsv --stat mean
```

Without `--stat`, per-base arrays are saved to a `.npz` file instead.

### Many TF folders at once

For a layout of `<input_root>/<TF>/6mer_sequences.tsv`:

```python
rs.build_datasets("human", "/path/to/human.bw",
                  input_root="data/", output_root="out/")
```

```bash
remap-signal build --species human --bigwig /path/to/human.bw \
                   --input-root data/ --output-root out/ \
                   [--tfs CTCF TP53] [--input-filename 6mer_sequences.tsv] [--seed 42]
```

Each TF gets its own output folder with the files from section 6.

## 8. Troubleshooting

| problem | cause and fix |
|---|---|
| `externally-managed-environment` | you used plain `pip`. Use `python -m pip` inside an activated environment. |
| install asks for a GitHub login | the repository is private. Make it public, or use a token. |
| `does not appear to be a Python project` | `pyproject.toml` is in a subfolder. Add `#subdirectory=<folder>` to the install URL. |
| `pyBigWig` fails to build | install it with conda: `conda install -c bioconda -c conda-forge pybigwig` |
| `No ReMap bigWig for human` | set `REMAP_HUMAN_BW`, or pass `bigwig="..."` |
| signal is all zeros | the coordinates and the bigWig use different assemblies, or the file has the wrong species. Check hg38 for human and mm10 for mouse. |
| warning about sequence length | `end - start` differs from the sequence length. The signal is cut or zero-padded to the sequence length. |

## 9. Develop

```bash
git clone https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git
cd Package_REMAP_SIGNAL
python -m pip install -e ".[dev]"
pytest
```

## 10. License

MIT
