"""Tests for bj_reference_builder.py (red/green TDD).

Builds the human (GRCh38) reference files BJ-Expression expects from public
sources: Ensembl 104 GTF, the NCBI no-alt analysis-set genome and Ensembl 104
cDNA/ncRNA FASTAs. Formats were read from the pipeline's own scripts:

- gtf_file: contigs must match the BAMs (chr-style, like the genome FASTA)
- tx2gene: 5 columns, no header, whitespace separated:
  transcript_id, transcript_biotype, gene_id, gene_biotype, gene_symbol
- genebody_ref: 12-column BED (geneBodyCoverageIsai.py skips anything else)
- bed_files_dir: <biotype>_transcripts.bed files (insert-size step)
"""

import gzip
import json
import importlib.util
import re
import sys
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
MODULE = SKILL_DIR / "bj_reference_builder.py"


def load():
    spec = importlib.util.spec_from_file_location("bj_reference_builder", MODULE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bj_reference_builder"] = mod
    spec.loader.exec_module(mod)
    return mod


def gtf_line(chrom, feature, start, end, strand="+", **attrs):
    a = " ".join(f'{k} "{v}";' for k, v in attrs.items())
    return "\t".join([chrom, "ensembl", feature, str(start), str(end), ".", strand, ".", a])


def tx_block(chrom="1", gene="ENSG1", tx="ENST1", exons=((100, 200), (300, 400)),
             strand="+", symbol="GENE1", gbio="protein_coding", tbio="protein_coding",
             gver="1", tver="2"):
    common = dict(gene_id=gene, gene_version=gver, transcript_id=tx,
                  transcript_version=tver, gene_name=symbol,
                  gene_biotype=gbio, transcript_biotype=tbio)
    lines = [gtf_line(chrom, "gene", exons[0][0], exons[-1][1], strand,
                      gene_id=gene, gene_version=gver, gene_name=symbol, gene_biotype=gbio),
             gtf_line(chrom, "transcript", exons[0][0], exons[-1][1], strand, **common)]
    for i, (s, e) in enumerate(exons, 1):
        lines.append(gtf_line(chrom, "exon", s, e, strand, exon_number=i, **common))
    return lines


# ──────────────────────────────────────────────
# Contig renaming
# ──────────────────────────────────────────────

@pytest.mark.parametrize("src,dst", [("1", "chr1"), ("22", "chr22"), ("X", "chrX"),
                                     ("Y", "chrY"), ("MT", "chrM")])
def test_primary_contigs_get_ucsc_names(src, dst):
    assert load().ucsc_name(src) == dst


@pytest.mark.parametrize("src", ["KI270728.1", "GL000194.1", "HSCHR1_CTG1_UNLOCALIZED", "chr1", ""])
def test_non_primary_contigs_have_no_name(src):
    assert load().ucsc_name(src) is None


def test_rename_gtf_renames_keeps_comments_and_counts_drops():
    mod = load()
    lines = ["#!genome-build GRCh38.p13", *tx_block("1"), *tx_block("MT", gene="ENSG2", tx="ENST2"),
             *tx_block("KI270728.1", gene="ENSG3", tx="ENST3")]
    out, stats = mod.rename_gtf(lines)
    assert out[0].startswith("#!")
    chroms = {ln.split("\t")[0] for ln in out if not ln.startswith("#")}
    assert chroms == {"chr1", "chrM"}
    assert stats["dropped_by_contig"] == {"KI270728.1": 4}
    assert stats["kept"] == 8


def test_rename_gtf_drops_contigs_missing_from_the_fasta():
    mod = load()
    lines = [*tx_block("1"), *tx_block("2", gene="ENSG2", tx="ENST2")]
    out, stats = mod.rename_gtf(lines, fasta_contigs={"chr1"})
    assert {ln.split("\t")[0] for ln in out} == {"chr1"}
    assert stats["not_in_fasta"] == {"chr2": 4}


@settings(deadline=None, max_examples=60)
@given(contigs=st.lists(st.sampled_from(["1", "2", "X", "MT", "KI270728.1", "GL000194.1"]),
                        min_size=0, max_size=12))
def test_property_every_data_line_is_kept_or_counted_as_dropped(contigs):
    mod = load()
    lines = []
    for i, c in enumerate(contigs):
        lines += tx_block(c, gene=f"ENSG{i}", tx=f"ENST{i}")
    out, stats = mod.rename_gtf(lines)
    dropped = sum(stats["dropped_by_contig"].values())
    assert len(out) == stats["kept"]
    assert stats["kept"] + dropped == len(lines)
    assert all(ln.split("\t")[0].startswith("chr") for ln in out)


# ──────────────────────────────────────────────
# Transcript records, tx2gene, BED12
# ──────────────────────────────────────────────

def test_transcript_record_has_exons_and_versions():
    mod = load()
    (rec,) = mod.transcripts(tx_block("chr1"))
    assert rec["tx_id"] == "ENST1.2" and rec["gene_id"] == "ENSG1.1"
    assert rec["exons"] == [(100, 200), (300, 400)]
    assert rec["symbol"] == "GENE1" and rec["tx_biotype"] == "protein_coding"


def test_versions_can_be_left_off():
    mod = load()
    (rec,) = mod.transcripts(tx_block("chr1"), with_versions=False)
    assert rec["tx_id"] == "ENST1" and rec["gene_id"] == "ENSG1"


def test_transcript_without_exons_is_skipped():
    mod = load()
    lines = [x for x in tx_block("chr1") if "\texon\t" not in x]
    assert mod.transcripts(lines) == []


def test_tx2gene_has_five_columns_in_the_order_the_pipeline_reads():
    mod = load()
    recs = mod.transcripts(tx_block("chr1", symbol="MT-ND1", tbio="protein_coding", gbio="protein_coding"))
    (line,) = mod.tx2gene_lines(recs)
    assert line.split("\t") == ["ENST1.2", "protein_coding", "ENSG1.1", "protein_coding", "MT-ND1"]


def test_missing_symbol_falls_back_to_the_gene_id():
    mod = load()
    lines = [re.sub(r' gene_name "[^"]*";', "", x) for x in tx_block("chr1")]
    (line,) = mod.tx2gene_lines(mod.transcripts(lines))
    assert line.split("\t")[4] == "ENSG1.1"


def test_bed12_matches_gtf_coordinates():
    mod = load()
    (rec,) = mod.transcripts(tx_block("chr1", exons=((101, 200), (301, 400)), strand="-"))
    f = mod.bed12(rec).split("\t")
    assert len(f) == 12
    assert f[:3] == ["chr1", "100", "400"]          # BED is 0-based half-open
    assert f[5] == "-" and f[9] == "2"
    assert f[10] == "100,100" and f[11] == "0,200"


@settings(deadline=None, max_examples=80)
@given(starts=st.lists(st.integers(1, 5000), min_size=1, max_size=8, unique=True),
       lens=st.lists(st.integers(1, 50), min_size=8, max_size=8),
       symbol=st.text(alphabet="ABCxyz09 -_.", min_size=1, max_size=12))
def test_property_bed12_blocks_and_tx2gene_fields_are_well_formed(starts, lens, symbol):
    mod = load()
    # non-overlapping exons: spread starts apart
    starts = sorted(s * 100 for s in starts)
    exons = tuple((s, s + lens[i] - 1) for i, s in enumerate(starts))
    recs = mod.transcripts(tx_block("chr1", exons=exons, symbol=symbol))
    assert len(recs) == 1
    f = mod.bed12(recs[0]).split("\t")
    assert len(f) == 12
    sizes = [int(x) for x in f[10].split(",")]
    offsets = [int(x) for x in f[11].split(",")]
    assert sizes == [e - s + 1 for s, e in exons]
    assert offsets[0] == 0 and offsets == sorted(offsets)
    assert int(f[1]) + offsets[-1] + sizes[-1] == int(f[2])
    assert int(f[9]) == len(exons)
    cols = mod.tx2gene_lines(recs)[0].split("\t")
    assert len(cols) == 5
    assert all(not re.search(r"\s", c) and c for c in cols), "read.table splits on whitespace"


# ──────────────────────────────────────────────
# Biotype BED files
# ──────────────────────────────────────────────

def test_one_bed_per_transcript_biotype_named_as_the_insert_size_step_reads(tmp_path):
    mod = load()
    recs = mod.transcripts([*tx_block("chr1", tx="ENST1", tbio="protein_coding"),
                            *tx_block("chr1", gene="ENSG2", tx="ENST2", tbio="lncRNA"),
                            *tx_block("chr1", gene="ENSG3", tx="ENST3", tbio="protein_coding")])
    written = mod.write_biotype_beds(recs, tmp_path / "biotypes")
    assert sorted(p.name for p in (tmp_path / "biotypes").iterdir()) == [
        "lncRNA_transcripts.bed", "protein_coding_transcripts.bed"]
    assert written == {"protein_coding": 2, "lncRNA": 1}
    lines = (tmp_path / "biotypes" / "protein_coding_transcripts.bed").read_text().splitlines()
    assert len(lines) == 2 and all(len(x.split("\t")) == 12 for x in lines)


def test_unsafe_biotype_names_are_refused_not_written(tmp_path):
    mod = load()
    recs = mod.transcripts(tx_block("chr1", tbio="bad/name"))
    with pytest.raises(mod.BuilderError, match="biotype"):
        mod.write_biotype_beds(recs, tmp_path / "biotypes")


def test_output_is_independent_of_input_order(tmp_path):
    mod = load()
    a = tx_block("chr1", tx="ENST1") + tx_block("chr1", gene="ENSG2", tx="ENST2", exons=((10, 20),))
    b = tx_block("chr1", gene="ENSG2", tx="ENST2", exons=((10, 20),)) + tx_block("chr1", tx="ENST1")
    assert mod.tx2gene_lines(mod.transcripts(a)) == mod.tx2gene_lines(mod.transcripts(b))
    assert [mod.bed12(r) for r in mod.transcripts(a)] == [mod.bed12(r) for r in mod.transcripts(b)]


# ──────────────────────────────────────────────
# FASTA helpers: transcripts and decoys
# ──────────────────────────────────────────────

def write_fa(path: Path, records: dict[str, str], gz=False):
    text = "".join(f">{h}\nACGT\n" for h in records.values()) if records else ""
    (gzip.open(path, "wt") if gz else open(path, "w")).write(text)


def test_fasta_ids_take_the_first_token(tmp_path):
    mod = load()
    p = tmp_path / "t.fa.gz"
    write_fa(p, {1: "ENST00000632684.1 cdna chromosome:GRCh38:7:1:2:1 gene:ENSG1.1", 2: "ENST2.3 cdna"}, gz=True)
    assert list(mod.fasta_ids(p)) == ["ENST00000632684.1", "ENST2.3"]


def test_transcript_fasta_is_filtered_to_ids_in_tx2gene(tmp_path):
    mod = load()
    a, b, out = tmp_path / "a.fa", tmp_path / "b.fa.gz", tmp_path / "out.fa"
    write_fa(a, {1: "ENST1.1 cdna", 2: "ENST2.1 cdna"})
    write_fa(b, {1: "ENST3.1 ncrna", 2: "ENST1.1 ncrna"}, gz=True)
    stats = mod.filter_fasta([a, b], {"ENST1.1", "ENST3.1"}, out)
    assert [h for h in mod.fasta_ids(out)] == ["ENST1.1", "ENST3.1"], "duplicates kept once"
    assert stats == {"kept": 2, "dropped": 1, "duplicates": 1}


def test_decoys_are_the_genome_contig_names(tmp_path):
    mod = load()
    p = tmp_path / "g.fa.gz"
    write_fa(p, {1: "chr1  AC:CM000663.2  LN:248956422", 2: "chrUn_KN707606v1_decoy  AC:KN707606.1"}, gz=True)
    assert mod.decoy_names(p) == ["chr1", "chrUn_KN707606v1_decoy"]


def test_tx2gene_coverage_check_reports_missing_ids():
    mod = load()
    assert mod.missing_from_tx2gene(["A", "B", "C"], {"A", "C"}) == ["B"]


# ──────────────────────────────────────────────
# Index commands (built, not run here)
# ──────────────────────────────────────────────

def test_star_command_matches_the_aligner_version_and_uses_the_renamed_gtf():
    mod = load()
    cmd = mod.star_index_command(Path("/w/genome.fa"), Path("/w/x.gtf"), Path("/w/idx"),
                                 threads=8, sjdb_overhang=100)
    s = " ".join(cmd)
    assert "--runMode genomeGenerate" in s and "--genomeFastaFiles /w/genome.fa" in s
    assert "--sjdbGTFfile /w/x.gtf" in s and "--sjdbOverhang 100" in s
    assert "--genomeDir /w/idx" in s and "--runThreadN 8" in s


def test_salmon_command_uses_decoys_and_the_gentrome():
    mod = load()
    cmd = mod.salmon_index_command(Path("/w/gentrome.fa"), Path("/w/decoys.txt"), Path("/w/idx"), threads=8)
    s = " ".join(cmd)
    assert s.startswith("salmon index") and "-t /w/gentrome.fa" in s
    assert "-d /w/decoys.txt" in s and "-i /w/idx" in s and "-p 8" in s


def test_verify_md5_detects_corruption(tmp_path):
    mod = load()
    f = tmp_path / "x.gz"
    f.write_bytes(b"hello")
    good = "5d41402abc4b2a76b9719d911017c592"
    assert mod.md5_matches(f, good) is True
    assert mod.md5_matches(f, "0" * 32) is False


# ──────────────────────────────────────────────
# End to end on tiny synthetic downloads
# ──────────────────────────────────────────────

import hashlib as _hashlib


def make_downloads(tmp_path: Path) -> Path:
    d = tmp_path / "downloads"
    d.mkdir()
    genome = d / "genome.fna.gz"
    with gzip.open(genome, "wt") as fh:
        for name in ("chr1", "chrM", "chrUn_X_decoy"):
            fh.write(f">{name}  AC:X  LN:100\nACGTACGT\n")
    md5 = _hashlib.md5(genome.read_bytes()).hexdigest()
    (d / "ncbi_md5checksums.txt").write_text(
        f"{md5}  ./GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz\n")
    gtf = [*tx_block("1", tx="ENST1", gene="ENSG1"),
           *tx_block("MT", tx="ENST2", gene="ENSG2", symbol="MT-ND1"),
           *tx_block("KI270728.1", tx="ENST3", gene="ENSG3")]
    with gzip.open(d / "ensembl104.gtf.gz", "wt") as fh:
        fh.write("#!genome-build GRCh38.p13\n" + "\n".join(gtf) + "\n")
    with gzip.open(d / "ensembl104.cdna.fa.gz", "wt") as fh:
        for tid in ("ENST1.2", "ENST3.2", "ENST9.1"):
            fh.write(f">{tid} cdna chromosome:GRCh38:1:1:2:1\nACGT\n")
    with gzip.open(d / "ensembl104.ncrna.fa.gz", "wt") as fh:
        fh.write(">ENST2.2 ncrna chromosome:GRCh38:MT:1:2:1\nACGT\n")
    return d


def test_prepare_writes_the_layout_the_pipeline_reads(tmp_path):
    mod = load()
    out = tmp_path / "refs"
    report = mod.prepare(make_downloads(tmp_path), out)
    seq = out / "genomes/Homo_sapiens/NCBI/GRCh38/Sequence"
    gtf = (seq / "star_chr_eqnames_wgs/Homo_sapiens.GRCh38.104.chr.gtf").read_text()
    assert "KI270728" not in gtf and "\nchr1\t" in "\n" + gtf and "\nchrM\t" in "\n" + gtf
    rows = (seq / "salmon/Ensemble_GRCh38_cDNA-ncRNA.complete.tx2gene").read_text().splitlines()
    assert [r.split("\t")[0] for r in rows] == ["ENST1.2", "ENST2.2"]
    assert all(len(r.split("\t")) == 5 for r in rows)
    assert sorted(p.name for p in (seq / "star_chr_eqnames_wgs/biotypes").iterdir()) == [
        "protein_coding_transcripts.bed"]
    work = out / "_work"
    assert list(mod.fasta_ids(work / "transcripts.fa")) == ["ENST1.2", "ENST2.2"], \
        "ENST3 (scaffold) and ENST9 (no annotation) must not reach Salmon"
    assert (work / "decoys.txt").read_text().split() == ["chr1", "chrM", "chrUn_X_decoy"]
    gent = list(mod.fasta_ids(work / "gentrome.fa"))
    assert gent == ["ENST1.2", "ENST2.2", "chr1", "chrM", "chrUn_X_decoy"], "transcripts first, then genome"
    assert report["transcript_ids_have_versions"] is True
    assert report["gtf"]["dropped_by_contig"] == {"KI270728.1": 4}


def test_every_salmon_transcript_is_in_tx2gene(tmp_path):
    """tximport fails on IDs missing from tx2gene; the filter must guarantee none."""
    mod = load()
    out = tmp_path / "refs"
    mod.prepare(make_downloads(tmp_path), out)
    rows = (out / mod.TX2GENE_PATH).read_text().splitlines()
    known = {r.split("\t")[0] for r in rows}
    assert mod.missing_from_tx2gene(mod.fasta_ids(out / "_work/transcripts.fa"), known) == []


def test_prepare_rejects_a_corrupt_genome_download(tmp_path):
    mod = load()
    d = make_downloads(tmp_path)
    (d / "genome.fna.gz").write_bytes((d / "genome.fna.gz").read_bytes() + b"x")
    with pytest.raises(mod.BuilderError, match="md5"):
        mod.prepare(d, tmp_path / "refs")


def test_prepare_names_missing_downloads(tmp_path):
    mod = load()
    d = make_downloads(tmp_path)
    (d / "ensembl104.ncrna.fa.gz").unlink()
    with pytest.raises(mod.BuilderError, match="ncrna"):
        mod.prepare(d, tmp_path / "refs")


def test_cli_prepare_runs_and_reports(tmp_path, capsys):
    mod = load()
    code = mod.main(["prepare", "--downloads", str(make_downloads(tmp_path)),
                     "--out", str(tmp_path / "refs")])
    assert code == 0 and "tx2gene_rows" in capsys.readouterr().out
    assert (tmp_path / "refs" / "prepare_report.json").exists()


# ──────────────────────────────────────────────
# Housekeeping-gene BED (genebody_ref)
# ──────────────────────────────────────────────
# The list carries gene versions from a different annotation than Ensembl 104
# (e.g. ENSG00000075624.17), so genes are matched on the UNVERSIONED gene ID.

def test_gene_list_strips_versions_and_ignores_blanks_and_comments():
    mod = load()
    text = "# housekeeping\nENSG1.17\tGENE1\n\n  ENSG2.5   GENE2  \nENSG1.9\tGENE1\n"
    assert mod.parse_gene_list(text) == [("ENSG1", "GENE1"), ("ENSG2", "GENE2")]


@pytest.mark.parametrize("bad", ["ENSG1.17\n", "ENSG1.17 A B\n", "notanid\tGENE\n"])
def test_gene_list_rejects_malformed_lines_with_a_line_number(bad):
    mod = load()
    with pytest.raises(mod.BuilderError, match="line 1"):
        mod.parse_gene_list(bad)


def _recs():
    mod = load()
    return mod.transcripts([
        *tx_block("chr1", gene="ENSG1", tx="ENST1", symbol="GENE1", gver="3"),
        *tx_block("chr1", gene="ENSG1", tx="ENST1b", symbol="GENE1", gver="3", exons=((1000, 1100),)),
        *tx_block("chr2", gene="ENSG2", tx="ENST2", symbol="GENE2", gver="1"),
        *tx_block("chr3", gene="ENSG3", tx="ENST3", symbol="GENE3", gver="1")])


def test_housekeeping_matches_on_unversioned_id_and_keeps_every_transcript():
    mod = load()
    lines, rep = mod.housekeeping_bed(_recs(), [("ENSG1", "GENE1"), ("ENSG2", "GENE2")])
    assert len(lines) == 3 and all(len(x.split("\t")) == 12 for x in lines)
    assert rep["per_gene_transcripts"] == {"ENSG1": 2, "ENSG2": 1}
    assert rep["unmatched"] == [] and rep["symbol_mismatch"] == [] and rep["symbol_only"] == []
    assert {x.split("\t")[3] for x in lines} == {"ENST1.2", "ENST1b.2", "ENST2.2"}


def test_unmatched_and_symbol_disagreement_are_reported_not_hidden():
    mod = load()
    lines, rep = mod.housekeeping_bed(_recs(), [("ENSG1", "WRONGNAME"), ("ENSG9", "NOPE")])
    assert rep["unmatched"] == ["ENSG9"]
    assert rep["symbol_mismatch"] == [["ENSG1", "WRONGNAME", "GENE1"]]
    assert len(lines) == 2, "an ID match is kept even when the symbol differs"


def test_symbol_only_match_is_excluded_but_reported():
    mod = load()
    lines, rep = mod.housekeeping_bed(_recs(), [("ENSG1", "GENE1"), ("ENSG77", "GENE2")])
    assert rep["symbol_only"] == [["ENSG77", "GENE2", "ENSG2"]]
    assert all(x.split("\t")[0] != "chr2" for x in lines)


def test_no_matches_is_an_error_not_an_empty_bed():
    mod = load()
    with pytest.raises(mod.BuilderError, match="no housekeeping genes"):
        mod.housekeeping_bed(_recs(), [("ENSG9", "NOPE")])


def test_housekeeping_command_writes_the_file_the_pipeline_names(tmp_path):
    mod = load()
    out = tmp_path / "refs"
    mod.prepare(make_downloads(tmp_path), out)
    genes = tmp_path / "genes.tsv"
    genes.write_text("ENSG1.17\tGENE1\nENSG2.4\tMT-ND1\n")
    assert mod.main(["housekeeping", "--genes", str(genes), "--out", str(out)]) == 0
    bed = out / "genomes/Homo_sapiens/NCBI/GRCh38/Sequence/star_chr_eqnames_wgs/hg38.HouseKeepingGenes.bed"
    rows = bed.read_text().splitlines()
    assert [r.split("\t")[3] for r in rows] == ["ENST1.2", "ENST2.2"]
    assert all(len(r.split("\t")) == 12 for r in rows)
    rep = json.loads((out / "housekeeping_report.json").read_text())
    assert rep["unmatched"] == []


def test_housekeeping_command_before_prepare_is_a_clear_error(tmp_path, capsys):
    mod = load()
    g = tmp_path / "g.tsv"
    g.write_text("ENSG1.1\tA\n")
    assert mod.main(["housekeeping", "--genes", str(g), "--out", str(tmp_path / "empty")]) == 1
    assert "prepare" in capsys.readouterr().err


# ──────────────────────────────────────────────
# build-indexes: STAR and Salmon inside the sandboxes, with a RAM guard
# ──────────────────────────────────────────────
# A human STAR index can exceed this container's RAM; an OOM kill restarts the
# container and loses background work, so we refuse to start without headroom.

import os as _os
import stat as _stat

FAKE_APPTAINER_IDX = f"""#!{sys.executable}
import json, os, sys
a = sys.argv[1:]
with open(os.environ["FAKE_IDX_LOG"], "a") as fh:
    fh.write(json.dumps({{"args": a, "cwd": os.getcwd()}}) + "\\n")
if os.environ.get("FAKE_IDX_FAIL"):
    sys.exit(0)  # pretend success but produce nothing
if "STAR" in a:
    d = a[a.index("--genomeDir") + 1]
    os.makedirs(d, exist_ok=True)
    for n in ("Genome", "SA", "SAindex", "genomeParameters.txt", "chrName.txt"):
        open(os.path.join(d, n), "w").write("x")
if "salmon" in a:
    d = a[a.index("-i") + 1]
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "info.json"), "w").write("{{}}")
"""

MEM_OK = "MemTotal: 70000000 kB\nMemFree: 1000 kB\nMemAvailable: 60000000 kB\n"   # ~57 GiB
MEM_LOW = "MemTotal: 37000000 kB\nMemFree: 1000 kB\nMemAvailable: 23000000 kB\n"  # ~22 GiB


def idx_env(tmp_path, monkeypatch):
    b = tmp_path / "bin"
    b.mkdir()
    exe = b / "apptainer"
    exe.write_text(FAKE_APPTAINER_IDX)
    exe.chmod(exe.stat().st_mode | _stat.S_IEXEC)
    log = tmp_path / "idx.log"
    monkeypatch.setenv("PATH", str(b))
    monkeypatch.setenv("FAKE_IDX_LOG", str(log))
    monkeypatch.delenv("FAKE_IDX_FAIL", raising=False)
    sb = tmp_path / "sandboxes"
    for im in ("quay.io-biocontainers-star-2.7.6a--0.img",
               "quay.io-biocontainers-salmon-1.6.0--h84f40af_0.img"):
        (sb / im).mkdir(parents=True)
    return log, sb


def prepared_refs(mod, tmp_path):
    out = tmp_path / "refs"
    mod.prepare(make_downloads(tmp_path), out)
    return out


def idx_calls(log):
    return [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []


def test_available_memory_is_parsed_from_meminfo():
    mod = load()
    assert abs(mod.available_memory_gb(MEM_OK) - 60000000 / 1024 / 1024) < 1e-6
    with pytest.raises(mod.BuilderError, match="MemAvailable"):
        mod.available_memory_gb("MemTotal: 1 kB\n")


@settings(deadline=None, max_examples=50)
@given(kb=st.integers(0, 10**9))
def test_property_memory_parse_roundtrips(kb):
    mod = load()
    assert mod.available_memory_gb(f"MemAvailable:   {kb} kB\n") == kb / 1024 / 1024


def test_low_memory_refuses_before_running_anything_and_mentions_the_resize(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    with pytest.raises(mod.BuilderError, match="resize"):
        mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_LOW)
    assert idx_calls(log) == [], "must not start STAR when it could be OOM-killed"
    assert not (out / mod.STAR_DIR / "index_ensemble").exists()


def test_both_indexes_are_built_in_their_sandboxes_with_the_prepared_files(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    rep = mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_OK)
    calls = idx_calls(log)
    star = next(c for c in calls if "STAR" in c["args"])["args"]
    salmon = next(c for c in calls if "salmon" in c["args"])["args"]
    assert star[:2] == ["exec", "--writable"] and "star-2.7.6a" in " ".join(star)
    assert str(out / mod.STAR_DIR / mod.GTF_NAME) in star and str(out / "_work/genome.fa") in star
    assert "salmon-1.6.0" in " ".join(salmon)
    assert str(out / "_work/gentrome.fa") in salmon and str(out / "_work/decoys.txt") in salmon
    assert (out / mod.STAR_DIR / "index_ensemble" / "Genome").is_file()
    assert (out / mod.SALMON_DIR / "index" / "info.json").is_file()
    assert rep["star"] == "built" and rep["salmon"] == "built"
    assert all(str(out) in " ".join(c["args"]) for c in calls)  # home/out is bound


def test_star_runs_in_the_work_dir_so_its_scratch_files_stay_under_the_output(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_OK)
    star_call = next(c for c in idx_calls(log) if "STAR" in c["args"])
    assert Path(star_call["cwd"]).resolve() == (out / "_work").resolve()


def test_existing_complete_indexes_are_skipped(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_OK)
    n = len(idx_calls(log))
    rep = mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_OK)
    assert rep["star"] == "skipped" and rep["salmon"] == "skipped"
    assert len(idx_calls(log)) == n


def test_an_index_that_comes_out_incomplete_is_an_error_not_a_success(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    monkeypatch.setenv("FAKE_IDX_FAIL", "1")
    with pytest.raises(mod.BuilderError, match="incomplete"):
        mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_OK)


def test_only_flag_builds_just_that_index(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    rep = mod.build_indexes(out, sb, threads=4, only="salmon", meminfo_text=MEM_LOW)
    assert rep["salmon"] == "built" and "star" not in rep, "Salmon needs less RAM than STAR"


def test_build_indexes_needs_prepare_first(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    with pytest.raises(mod.BuilderError, match="prepare"):
        mod.build_indexes(tmp_path / "empty", sb, threads=4, meminfo_text=MEM_OK)


def test_missing_sandbox_names_the_image_and_the_fix(tmp_path, monkeypatch):
    mod = load()
    log, sb = idx_env(tmp_path, monkeypatch)
    out = prepared_refs(mod, tmp_path)
    import shutil as _sh
    _sh.rmtree(sb / "quay.io-biocontainers-star-2.7.6a--0.img")
    with pytest.raises(mod.BuilderError, match="apptainer_sandboxes"):
        mod.build_indexes(out, sb, threads=4, meminfo_text=MEM_OK)


# ──────────────────────────────────────────────
# Cell-typing reference (celltype_ref)
# ──────────────────────────────────────────────
# The pipeline's script runs readRDS("db_celltype_qbio.RDS") on the staged file,
# so whatever is wired in must be named exactly that, at
# <genomes_base>/dev-resources/rnaseq_celltype_qbio/. It also indexes the list by
# five fixed entry names, so other lists "work" only until that step fails.

def test_the_five_entry_names_the_script_hardcodes():
    mod = load()
    assert mod.CELLTYPE_EXPECTED == ("HPCA_singleR", "GTEx_singleR", "GTExmeds_singleR",
                                     "TCGA_STN_singleR", "TCGA_PT_singleR")


@settings(deadline=None, max_examples=60)
@given(names=st.lists(st.sampled_from(["HPCA_singleR", "GTEx_singleR", "GTExmeds_singleR",
                                       "TCGA_STN_singleR", "TCGA_PT_singleR", "HPCA", "other"]),
                      max_size=7))
def test_property_missing_entries_are_exactly_the_expected_ones_not_present(names):
    mod = load()
    missing = mod.missing_celltype_entries(names)
    assert set(missing) == set(mod.CELLTYPE_EXPECTED) - set(names)
    assert list(missing) == [n for n in mod.CELLTYPE_EXPECTED if n in missing], "stable order"


def test_celltype_ref_is_copied_under_the_exact_name_the_script_reads(tmp_path):
    mod = load()
    src = tmp_path / "celldex_HPCA_plus_GTEx.rds"
    src.write_bytes(b"rds-bytes")
    rep = mod.place_celltype_ref(tmp_path / "refs", src,
                                 names=["HPCA", "GTEx_v9_tissue_aggregate_min10"])
    dest = tmp_path / "refs/dev-resources/rnaseq_celltype_qbio/db_celltype_qbio.RDS"
    assert dest.read_bytes() == b"rds-bytes" and not dest.is_symlink()
    assert rep["sha256"] == __import__("hashlib").sha256(b"rds-bytes").hexdigest()
    assert rep["source"] == str(src) and rep["placed"] == str(dest)
    assert rep["compatible_with_pipeline_script"] is False
    assert rep["expected_entries_missing"] == list(mod.CELLTYPE_EXPECTED)


def test_a_list_with_all_five_entries_is_reported_compatible(tmp_path):
    mod = load()
    src = tmp_path / "x.rds"
    src.write_bytes(b"1")
    rep = mod.place_celltype_ref(tmp_path / "refs", src, names=list(mod.CELLTYPE_EXPECTED))
    assert rep["compatible_with_pipeline_script"] is True and rep["expected_entries_missing"] == []


def test_unknown_entry_names_are_reported_as_unchecked_not_compatible(tmp_path):
    mod = load()
    src = tmp_path / "x.rds"
    src.write_bytes(b"1")
    rep = mod.place_celltype_ref(tmp_path / "refs", src, names=None)
    assert rep["compatible_with_pipeline_script"] is None and rep["inspected"] is False


def test_a_different_existing_file_is_not_overwritten_without_force(tmp_path):
    mod = load()
    a, b = tmp_path / "a.rds", tmp_path / "b.rds"
    a.write_bytes(b"AAA")
    b.write_bytes(b"BBB")
    mod.place_celltype_ref(tmp_path / "refs", a)
    mod.place_celltype_ref(tmp_path / "refs", a)  # same content: fine, idempotent
    with pytest.raises(mod.BuilderError, match="force"):
        mod.place_celltype_ref(tmp_path / "refs", b)
    dest = tmp_path / "refs/dev-resources/rnaseq_celltype_qbio/db_celltype_qbio.RDS"
    assert dest.read_bytes() == b"AAA"
    mod.place_celltype_ref(tmp_path / "refs", b, force=True)
    assert dest.read_bytes() == b"BBB"


def test_missing_celltype_source_is_named(tmp_path):
    mod = load()
    with pytest.raises(mod.BuilderError, match="nope.rds"):
        mod.place_celltype_ref(tmp_path / "refs", tmp_path / "nope.rds")


def test_celltype_command_writes_a_report_and_says_when_it_could_not_inspect(tmp_path, capsys, monkeypatch):
    mod = load()
    monkeypatch.setenv("PATH", str(tmp_path))  # no apptainer: inspection must be skipped, not crash
    src = tmp_path / "x.rds"
    src.write_bytes(b"1")
    code = mod.main(["celltype", "--rds", str(src), "--out", str(tmp_path / "refs")])
    assert code == 0
    rep = json.loads((tmp_path / "refs" / "celltype_report.json").read_text())
    assert rep["inspected"] is False and rep["compatible_with_pipeline_script"] is None
    assert "not inspected" in capsys.readouterr().out.lower()


def test_celltype_command_inspects_the_placed_copy_not_the_source(tmp_path, monkeypatch, capsys):
    """The source may be on a mount the container cannot see (/orgbucket); the
    copy under the references root is the one the container can read."""
    mod = load()
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "inspect.log"
    exe = bindir / "apptainer"
    exe.write_text(f"#!{sys.executable}\nimport sys\nopen({str(log)!r}, 'a').write(' '.join(sys.argv[1:]) + '\\n')\n"
                   "print('HPCA'); print('GTEx_v9_tissue_aggregate_min10')\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir))
    sandboxes = tmp_path / "sandboxes"
    (sandboxes / mod._cache_name(mod.CELLTYPE_IMAGE)).mkdir(parents=True)
    src = tmp_path / "elsewhere" / "ref.rds"
    src.parent.mkdir()
    src.write_bytes(b"1")
    code = mod.main(["celltype", "--rds", str(src), "--out", str(tmp_path / "refs"),
                     "--sandboxes", str(sandboxes)])
    assert code == 0
    called = log.read_text()
    placed = str(tmp_path / "refs/dev-resources/rnaseq_celltype_qbio/db_celltype_qbio.RDS")
    assert placed in called and str(src) not in called
    rep = json.loads((tmp_path / "refs" / "celltype_report.json").read_text())
    assert rep["inspected"] is True and rep["compatible_with_pipeline_script"] is False
    assert rep["entries_found"] == ["HPCA", "GTEx_v9_tissue_aggregate_min10"]
    assert "WARNING" in capsys.readouterr().out
