# Package: REMAP_SIGNAL

### Instructions to get this package into your environment
#### 1. Activate your environment
#### 2. `conda install -c bioconda -c conda-forge pybigwig scikit-learn pandas numpy tqdm -y`
#### 3. `python -m pip install git+https://github.com/NPAPINENI/Package_REMAP_SIGNAL.git`

---

## What it does

Give it sequences with their genomic positions and a species. It adds the **ReMap signal**
for every base and splits the data into **train / dev / test**.

## What you give

A TSV or CSV file with these columns:

| chr | start | end | label | sequence |
|---|---|---|---|---|
| chr1 | 1000000 | 1000300 | 1 | ACGT... |
| chr2 | 5000000 | 5000300 | 0 | TTGA... |

(`start`/`end` are 0-based; `end - start` = sequence length.)

## What you get

In your output folder:

| file | content |
|---|---|
| `train.tsv`, `dev.tsv`, `test.tsv` | `Sequence`, `Label` |
| `sequence_train.tsv`, `sequence_dev.tsv`, `sequence_test.tsv` | ReMap signal per base (`nucleotide_1 ... nucleotide_N`), same row order as above |

## How to use

**Human and mouse are built in:**

```python
import remap_signal as rs

rs.make_splits("my_file.tsv", species="human", output_dir="out_human/")   # hg38
rs.make_splits("my_file.tsv", species="mouse", output_dir="out_mouse/")   # mm39
```

or from the terminal:

```bash
remap-signal split --input my_file.tsv --species human --output-dir out_human/
```

**Any other species:** give your own bigWig file with the signal.

```python
rs.make_splits("my_file.tsv", species="zebrafish", output_dir="out_fish/",
               bigwig="/path/to/zebrafish_signal.bw")
```

```bash
remap-signal split --input my_file.tsv --species zebrafish --output-dir out_fish/ \
                   --bigwig /path/to/zebrafish_signal.bw
```

The genome version of your coordinates must match the species (human hg38, mouse mm39) or
your bigWig.
