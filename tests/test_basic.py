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


def test_missing_explicit_bigwig_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        rs.resolve_bigwig("human", str(tmp_path / "nope.bw"))


def test_species_aliases():
    from remap_signal.core import _norm_species
    assert _norm_species("human") == "hg38" and _norm_species("mouse") == "mm39"
    assert _norm_species("mm10") == "mm10" and _norm_species("DM6") == "dm6"
    with pytest.raises(ValueError):
        _norm_species("zebrafish")


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


def _gz(tmp_path, text):
    import gzip
    f = tmp_path / "x.bed.gz"
    with gzip.open(f, "wt") as h:
        h.write(text)
    return f


def test_build_crm_interval_density(tmp_path, monkeypatch):
    bed = _gz(tmp_path, "chr1\t100\t200\tA\nchr1\t150\t250\tB\nchr2\t10\t20\tC\n")
    monkeypatch.setenv("REMAP_CACHE_DIR", str(tmp_path / "cache"))
    p = rs.download("mouse", kind="crm", url="file://" + str(bed))
    with rs.RemapSignal("mouse", p) as r:
        v = r.get_signal("chr1", 90, 260)[0]
        assert v[:10].sum() == 0 and (v[10:60] == 1).all()
        assert (v[60:110] == 2).all() and (v[110:160] == 1).all() and v[160:].sum() == 0
        assert r.get_signal("chr2", 10, 20)[0].tolist() == [1.0] * 10


def test_build_all_uses_summit_plus_minus_25(tmp_path, monkeypatch):
    # BED columns: chr start end name score strand thickStart thickEnd
    bed = _gz(tmp_path, "chr1\t0\t500\tA\t0\t.\t100\t101\nchr1\t0\t500\tB\t0\t.\t120\t121\n")
    monkeypatch.setenv("REMAP_CACHE_DIR", str(tmp_path / "cache"))
    p = rs.download("mm10", kind="all", url="file://" + str(bed))
    with rs.RemapSignal("mm10", p) as r:
        v = r.get_signal("chr1", 70, 150)[0]
        assert v[:5].sum() == 0                       # before 75
        assert (v[5:25] == 1).all()                   # 75..94: peak A only
        assert (v[25:56] == 2).all()                  # 95..125: both
        assert (v[56:76] == 1).all() and v[76:].sum() == 0   # 126..145: peak B only


def test_direct_download_validates_bigwig(tmp_path, monkeypatch, bw):
    import sys
    dl = sys.modules["remap_signal.download"]
    monkeypatch.setattr(dl, "DIRECT", {"hg38": {"url": "file://" + bw, "md5": None}})
    monkeypatch.setenv("REMAP_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.delenv("REMAP_HUMAN_BW", raising=False)
    monkeypatch.delenv("REMAP_HG38_BW", raising=False)
    p = rs.download("human")
    assert open(p, "rb").read() == open(bw, "rb").read()
    assert rs.resolve_bigwig("human") == p            # cached afterwards
    junk = tmp_path / "junk.bw"
    junk.write_text("<html>not a bigwig</html>")
    monkeypatch.setattr(dl, "DIRECT", {"hg38": {"url": "file://" + str(junk), "md5": None}})
    with pytest.raises(OSError):
        rs.download("human", force=True)
    assert open(p, "rb").read() == open(bw, "rb").read()   # good cache survives a bad re-download


def test_other_species_with_own_bigwig(bw, tmp_path):
    rows = [dict(chr="chr1", start=900 + i * 10, end=1000 + i * 10, label=i % 2,
                 sequence="ACGT" * 25) for i in range(40)]
    f = tmp_path / "fish.tsv"
    pd.DataFrame(rows).to_csv(f, sep="\t", index=False)
    res = rs.make_splits(f, "zebrafish", tmp_path / "o", bigwig=bw)
    assert len(res["train"]) == 28
    with pytest.raises(ValueError):                     # unknown species without a bigWig
        rs.make_splits(f, "zebrafish", tmp_path / "o2")


def test_friendly_helpers(bw):
    s = rs.get_signal("human", "chr1", 950, 1150, bigwig=bw)
    assert s.shape == (200,) and (s[50:150] == 4).all() and s[:50].sum() == 0
    df = pd.DataFrame(dict(chr=["chr1", "1"], start=[1000, 1050], end=[1010, 1070]))
    out = rs.add_signal(df, "zebrafish", bigwig=bw)
    assert out.shape == (2, 3 + 20) and (out["nucleotide_1"] == 4).all()
    assert out.loc[0, "nucleotide_11":].sum() == 0       # shorter row is zero-padded
