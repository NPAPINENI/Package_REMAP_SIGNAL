import numpy as np
import pandas as pd
import pyBigWig
import pytest

import remap_signal as rs


@pytest.fixture
def bw(tmp_path):
    p = str(tmp_path / "x.bw")
    f = pyBigWig.open(p, "w")
    f.addHeader([("chr1", 10000)])
    f.addEntries("chr1", [1000], values=[4.0], span=100, step=100)
    f.close()
    return p


def test_signal_values_and_alignment(bw):
    with rs.RemapSignal("human", bw) as r:
        a = r.get_signal("chr1", 950, 1150, "A" * 200)
    assert a.shape == (1, 200)
    assert a[0, :50].sum() == 0 and (a[0, 50:150] == 4).all() and a[0, 150:].sum() == 0


def test_chrom_naming_and_edges(bw):
    with rs.RemapSignal("mouse", bw) as r:
        assert r.get_signal("1", 1000, 1010).shape == (1, 10)
        assert r.get_signal("chrZ", 0, 10).sum() == 0
        assert r.get_signal("chr1", -5, 5).shape == (1, 10)


def test_length_mismatch_raises(bw):
    with rs.RemapSignal("human", bw) as r:
        with pytest.raises(ValueError):
            r.get_signal("chr1", 0, 10, "ACG")


def test_missing_bigwig_message(monkeypatch, tmp_path):
    monkeypatch.delenv("REMAP_HUMAN_BW", raising=False)
    f = tmp_path / "in.tsv"
    pd.DataFrame(dict(chr=["chr1"], start=[0], end=[4], label=[0], sequence=["ACGT"])).to_csv(f, sep="\t", index=False)
    with pytest.raises(ValueError, match="REMAP_HUMAN_BW"):
        rs.make_splits(f, "human", tmp_path / "o")


def test_make_splits(bw, tmp_path):
    rows = [dict(chr="chr1", start=900 + i * 10, end=1000 + i * 10, label=i % 2,
                 sequence="ACGT" * 25) for i in range(40)]
    f = tmp_path / "in.tsv"
    pd.DataFrame(rows).to_csv(f, sep="\t", index=False)
    out = tmp_path / "out"
    res = rs.make_splits(f, "human", out, bigwig=bw)
    assert (len(res["train"]), len(res["dev"]), len(res["test"])) == (28, 6, 6)
    for s in ("train", "dev", "test"):
        assert (out / f"{s}.tsv").exists()
        sig = pd.read_csv(out / f"sequence_{s}.tsv", sep="\t")
        assert sig.shape[1] == 100 and len(sig) == len(pd.read_csv(out / f"{s}.tsv", sep="\t"))
    assert np.isfinite(res["all"].iloc[:, 2:].values).all()
