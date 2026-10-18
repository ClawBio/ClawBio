"""Single-position reads from a bgzipped, tabix-indexed sumstats file (local or HTTPS).

Why this is not `pysam.TabixFile.fetch`:

1. **Pan-UKBB's indexes defeat a narrow htslib query.** Measured 2026-10-17 on
   `biomarkers-30780-both_sexes-irnt` (the LDL file): an exact-position
   `TabixFile.fetch(chrom, pos - 1, pos)` returned nothing for 16 of 250 randomly sampled
   positions that a +/-500 kb fetch of the same file shows are present (chr5:50-52 Mb),
   and for rs12740374 itself (chr1:109817590) in 32 of 32 randomly sampled EUR outcome
   files. FinnGen R12 (`I9_CHD`, 200 sampled positions) showed no misses. The cause, read
   from the index bytes of `icd10-I62-both_sexes`: the index's virtual offsets point at an
   EMPTY bgzip block (the end-of-partition marker the writer leaves between concatenated
   parts) with an in-block offset that only makes sense in the block AFTER it, so htslib
   seeks into a zero-length block and returns no records. A miss would have been served
   as "variant absent", the silent-absence defect. This reader walks the index's linear
   offsets itself and decompresses the byte stream as one concatenation, in which an empty
   block contributes zero bytes, so the same offset lands on the right record. Checked
   against a +/-500 kb htslib fetch as truth (2026-10-17): agreement on all 32 files at
   rs12740374, and on every sampled position, present and absent: 250 + 62 (LDL, chr5),
   150 + 37 (`icd10-I62`, chr1), 150 + 37 (FinnGen `I9_CHD`, chr1).
2. **htslib writes a downloaded remote index into the current directory.** An on-demand
   lookup over the 7,160 Pan-UKBB EUR outcomes would leave ~15.8 GB of `.tbi` files
   wherever the process was started (7,160 x the 2.2 MB index of the LDL file; an
   estimate, not a measured total). Indexes here go to an explicit cache directory.

The reader is format-agnostic: the header line is read from the file itself and columns are
selected by name downstream.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import struct
import tempfile
import threading
import zlib
from dataclasses import dataclass

import requests

TILE_SHIFT = 14  # tabix linear-index tile = 16 kb
DEFAULT_TIMEOUT_S = 60.0
_INITIAL_READ = 128 * 1024
_MAX_READ = 8 * 1024 * 1024
INDEX_CACHE_ENV = "SNP_PHEWAS_INDEX_CACHE"


RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
RETRY_BACKOFF_S = 1.0


def get_with_retry(session: "requests.Session", url: str, *, attempts: int = 3,
                   backoff_s: float | None = None, **kw) -> "requests.Response":
    """GET with a short bounded retry on a transport error or a 429/5xx answer.

    For the single-request calls a whole dataset depends on (the assembly mapper, the
    FinnGen variant API): one transient 503 there used to mark every outcome of the
    dataset unreachable (measured 2026-10-17: Ensembl answered 503 once for a
    589-outcome lookup that succeeded on the next try). The last answer or error is
    returned or re-raised unchanged, so the caller's classification is untouched."""
    import time

    for i in range(attempts):
        try:
            r = session.get(url, **kw)
        except requests.RequestException:
            if i == attempts - 1:
                raise
        else:
            if r.status_code not in RETRY_STATUSES or i == attempts - 1:
                return r
        time.sleep((RETRY_BACKOFF_S if backoff_s is None else backoff_s) * (2 ** i))
    raise AssertionError("unreachable")  # pragma: no cover


class SourceUnreachable(ConnectionError):
    """The source could not be READ (network, timeout, 5xx, throttling).

    Never means the variant is absent. A ConnectionError subclass so a transport-failure
    `isinstance` check recognises it, and its own type so it cannot be coerced into an
    empty result by a caller that only checks for emptiness.
    """


class SourceFileMissing(FileNotFoundError):
    """The source answered, and the file is not there (HTTP 403/404, or no local file)."""


class IndexFormatError(ValueError):
    """The index bytes are not a tabix index this reader understands."""


@dataclass(frozen=True)
class TabixMeta:
    names: tuple[str, ...]
    col_seq: int
    col_beg: int
    meta_char: str
    skip: int


_tls = threading.local()

# Index-cache accounting for the run report: how many remote indexes were downloaded (and
# how many bytes) versus found already cached.
_index_lock = threading.Lock()
INDEX_STATS = {"downloaded": 0, "bytes_downloaded": 0, "cache_hits": 0}


def reset_index_stats() -> None:
    with _index_lock:
        for k in INDEX_STATS:
            INDEX_STATS[k] = 0


def _count(key: str, n: int = 1) -> None:
    with _index_lock:
        INDEX_STATS[key] += n


def _session() -> requests.Session:
    s = getattr(_tls, "session", None)
    if s is None:
        s = requests.Session()
        _tls.session = s
    return s


def _is_url(path: str) -> bool:
    return path.startswith(("http://", "https://"))


def _http_get(url: str, *, byte_range: tuple[int, int] | None, timeout_s: float) -> bytes:
    headers = {"Range": f"bytes={byte_range[0]}-{byte_range[1]}"} if byte_range else {}
    try:
        r = _session().get(url, headers=headers, timeout=timeout_s)
    except requests.RequestException as e:
        raise SourceUnreachable(f"could not read {url}: {type(e).__name__}") from e
    if r.status_code in (403, 404):
        raise SourceFileMissing(f"{url} answered HTTP {r.status_code}")
    if r.status_code == 416:  # range past EOF: nothing more to read
        return b""
    if r.status_code >= 400:
        raise SourceUnreachable(f"{url} answered HTTP {r.status_code}")
    return r.content


def read_bytes(path: str, start: int, length: int, *, timeout_s: float = DEFAULT_TIMEOUT_S) -> bytes:
    """`length` bytes from `start` of a local file or an HTTPS object (fewer at EOF)."""
    if _is_url(path):
        return _http_get(path, byte_range=(start, start + length - 1), timeout_s=timeout_s)
    try:
        with open(path, "rb") as fh:
            fh.seek(start)
            return fh.read(length)
    except FileNotFoundError as e:
        raise SourceFileMissing(str(e)) from e


def _inflate_blocks(data: bytes) -> tuple[bytes, int]:
    """Inflate every COMPLETE bgzip block in `data`; return (payload, bytes consumed).

    Empty blocks inflate to b"" and so vanish from the concatenation, which is what makes
    an offset into "the block after the empty one" land correctly.
    """
    out = bytearray()
    i = 0
    n = len(data)
    while i + 18 <= n:
        if data[i:i + 4] != b"\x1f\x8b\x08\x04":
            raise IndexFormatError(f"not a bgzip block at byte {i}")
        bsize = struct.unpack_from("<H", data, i + 16)[0] + 1
        if i + bsize > n:
            break
        out += zlib.decompress(data[i + 18:i + bsize - 8], -15)
        i += bsize
    return bytes(out), i


def parse_index(raw: bytes) -> tuple[TabixMeta, dict[str, list[int]]]:
    """Parse a (gzip-compressed) tabix index into its meta + per-contig linear offsets."""
    b = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    if b[:4] != b"TBI\x01":
        raise IndexFormatError("missing TBI magic")
    n_ref, _fmt, col_seq, col_beg, _col_end, meta, skip, l_nm = struct.unpack_from("<8i", b, 4)
    o = 36
    names = tuple(x.decode() for x in b[o:o + l_nm].split(b"\0")[:-1])
    o += l_nm
    linear: dict[str, list[int]] = {}
    for name in names[:n_ref]:
        (n_bin,) = struct.unpack_from("<i", b, o)
        o += 4
        for _ in range(n_bin):
            _bin, n_chunk = struct.unpack_from("<Ii", b, o)
            o += 8 + 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", b, o)
        o += 4
        linear[name] = list(struct.unpack_from(f"<{n_intv}Q", b, o))
        o += 8 * n_intv
    return TabixMeta(names, col_seq, col_beg, chr(meta), skip), linear


def default_index_cache() -> str:
    return os.environ.get(INDEX_CACHE_ENV) or os.path.join(
        os.path.expanduser("~"), ".clawbio", "snp_phewas", "index"
    )


def ensure_index(index_path: str, cache_dir: str | None = None, *,
                 timeout_s: float = DEFAULT_TIMEOUT_S) -> str:
    """Local path of the index: a local index is used in place, a remote one is fetched
    once into `cache_dir` (atomic rename, so a half-written index is never read)."""
    if not _is_url(index_path):
        if not os.path.exists(index_path):
            raise SourceFileMissing(index_path)
        return index_path
    cache_dir = cache_dir or default_index_cache()
    os.makedirs(cache_dir, exist_ok=True)
    # Keyed by the full URL, not the basename: two sources (or two releases at different
    # paths) can publish the same file name. Callers pass a release-scoped `cache_dir`.
    key = hashlib.sha256(index_path.encode()).hexdigest()[:16]
    local = os.path.join(cache_dir, f"{key}-{index_path.rsplit('/', 1)[-1]}")
    if os.path.exists(local):
        _count("cache_hits")
        return local
    data = _http_get(index_path, byte_range=None, timeout_s=timeout_s)
    _count("downloaded")
    _count("bytes_downloaded", len(data))
    fd, tmp = tempfile.mkstemp(dir=cache_dir, suffix=".part")
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
    os.replace(tmp, local)
    return local


def read_header(data_path: str, *, timeout_s: float = DEFAULT_TIMEOUT_S) -> list[str]:
    """The file's column names: its first line, with a leading `#` stripped."""
    raw = read_bytes(data_path, 0, _INITIAL_READ, timeout_s=timeout_s)
    payload, _ = _inflate_blocks(raw)
    line = payload.split(b"\n", 1)[0].decode()
    return line.lstrip("#").split("\t")


def _contig(meta: TabixMeta, chrom: str) -> str | None:
    c = chrom[3:] if chrom.startswith("chr") else chrom
    alias = {"X": "23", "Y": "24", "MT": "25", "M": "25"}.get(c, "")
    for cand in (c, f"chr{c}", alias):
        if cand and cand in meta.names:
            return cand
    return None


def fetch_position(data_path: str, index_local: str, chrom: str, pos: int, *,
                   timeout_s: float = DEFAULT_TIMEOUT_S) -> list[list[str]]:
    """Every record at exactly (`chrom`, `pos`) (1-based), as split fields.

    Starts at the linear-index offset of the tile BEFORE the target's (one tile of
    slack against an off-by-one index) and scans forward until a record passes `pos`.
    An empty list means the file was read past the position and holds no record there.
    """
    with open(index_local, "rb") as fh:
        meta, linear = parse_index(fh.read())
    contig = _contig(meta, chrom)
    if contig is None:
        return []
    lin = linear.get(contig, [])
    tile = (pos - 1) >> TILE_SHIFT
    if tile >= len(lin):
        cands = [v for v in lin if v]
        if not cands:
            return []
        voff = cands[-1]
    else:
        cands = [v for v in lin[max(0, tile - 1):tile + 1] if v]
        if not cands:
            # Both tiles empty in the index: start from the nearest earlier tile that has an
            # offset, never byte 0 of a large file (a scan from 0 would hit the read cap and
            # mislabel a real present-or-absent answer as an outage).
            earlier = [v for v in lin[:max(0, tile - 1)] if v]
            cands = earlier[-1:]
        voff = min(cands) if cands else 0
    coff, uoff = voff >> 16, voff & 0xFFFF
    seq_i, beg_i = meta.col_seq - 1, meta.col_beg - 1

    want = _INITIAL_READ
    while True:
        raw = read_bytes(data_path, coff, want, timeout_s=timeout_s)
        payload, consumed = _inflate_blocks(raw)
        at_eof = len(raw) < want
        text = payload[uoff:]
        lines = text.split(b"\n")
        if not at_eof:
            lines = lines[:-1]  # the last line may be cut mid-record
        hits: list[list[str]] = []
        passed = False
        for ln in lines:
            if not ln or ln[:1] == meta.meta_char.encode():
                continue
            f = ln.decode().split("\t")
            if len(f) <= max(seq_i, beg_i) or not f[beg_i].isdigit():
                continue  # a header line (Pan-UKBB's `chr pos ...` is not `#`-prefixed)
            if f[seq_i] != contig:
                if hits or passed or f[seq_i] in meta.names[meta.names.index(contig) + 1:]:
                    passed = True
                    break
                continue
            p = int(f[beg_i])
            if p == pos:
                hits.append(f)
            elif p > pos:
                passed = True
                break
        if passed or at_eof:
            return hits
        if want >= _MAX_READ:
            raise SourceUnreachable(
                f"read {want} bytes of {data_path} from the index offset without reaching "
                f"{chrom}:{pos}; refusing to report the variant as absent"
            )
        want *= 2
