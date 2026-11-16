# REMAP_SIGNAL


Attach **ReMap** regulatory signal to genomic sequences and generate stratified
**train(70%)/ dev(15%) / test(15%)** datasets, in a single call.

- Built-in support for **human (hg38)** and **mouse (mm39)**
- Works with **any other species** through a user-supplied bigWig
- Per-base signal, aligned to each input sequence
- Python API and command-line interface

---

## Installation

```bash
# 1. activate your environment
conda activate <your-env>

# 2. install dependencies
conda install -c bioconda -c conda-forge pybigwig scikit-learn pandas numpy tqdm -y

# 3. install the package
python -m pip install git+https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git
```

Verify the installation:

```bash
remap-signal --help
```

**Requirements:** Python 3.8 or newer.

---

## Quick start

```python
import remap_signal as rs

rs.make_splits("my_file.tsv", species="human", output_dir="out_human/")
rs.make_splits("my_file.tsv", species="mouse", output_dir="out_mouse/")
```

Command line:

```bash
remap-signal split --input my_file.tsv --species human --output-dir out_human/
```

---

## Input

A tsv/csv file with one row per sequence:

| Column | Description |
|---|---|
| `chr` | Chromosome (`chr1` or `1`) |
| `start` | Start position, 0-based |
| `end` | End position, exclusive |
| `label` | Class label (used to stratify the split) |
| `sequence` | DNA sequence (`end - start` bases) |

Additional columns are ignored. Coordinates must match the genome version of the selected
species.

## Output

Written to `output_dir`:

| File | Content |
|---|---|
| `train.tsv`, `dev.tsv`, `test.tsv` | `Sequence`, `Label` |
| `sequence_train.tsv`, `sequence_dev.tsv`, `sequence_test.tsv` | ReMap signal per base (`nucleotide_1` … `nucleotide_N`); row *i* corresponds to row *i* of the matching split file |
| `Remap_signal_<name>.tsv` | Full table before splitting |

---

## Supported species

| `species=` | Genome |
|---|---|
| `"human"` | hg38 |
| `"mouse"` | mm39 |
| any other name | provide `bigwig=` (see below) |

### Other species or your own signal track

```python
rs.make_splits("my_file.tsv", species="zebrafish", output_dir="out_fish/",
               bigwig="/path/to/zebrafish_signal.bw")
```

```bash
remap-signal split --input my_file.tsv --species zebrafish --output-dir out_fish/ \
                   --bigwig /path/to/zebrafish_signal.bw
```

A bigWig passed with `bigwig=` is also used for human and mouse when provided.

---

## Data sources and citation

ReMap 2022 density tracks are provided by the ReMap project. If you use this package,
please cite:

> Hammal F, de Langen P, Bergon A, Lopez F, Ballester B. **ReMap 2022: a database of Human,
> Mouse, Drosophila and Arabidopsis regulatory regions from an integrative analysis of
> DNA-binding sequencing experiments.** *Nucleic Acids Research* 2022;50(D1):D316–D325.
> https://doi.org/10.1093/nar/gkab996

## Development

```bash
git clone https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git
cd Package_REMAP_SIGNAL
python -m pip install -e ".[dev]"
pytest
```

## License

MIT © Nimisha Papineni
