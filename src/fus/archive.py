"""Acoustic recordings read straight out of the lab's Dropbox archive, without unpacking it.

``data/acoustic/US_Data.zip`` (48.6 GB, from Dropbox 2026-10-19) holds all four controller
sessions, one folder per day (``061126_FUS_AAV_Controller`` = 2026-06-11, mice 1-6; then
06-16 mice 7-12, 07-01 mice 13-18, 07-07 mice 19-24). Members are *stored*, not
compressed, so one recording can be copied out in well under a second. **The zip is the only
copy**: nothing is unpacked (2026-10-19, to keep one copy of the data). Analyses go through
:func:`local_path` / :func:`paths`. They yield a temporary copy that is deleted afterwards,
or a file in ``data/acoustic/<yyyymmdd>/`` if someone has unpacked one there at the same size.

Only ``Mouse_Cntr_NN_Target*.mat`` are listed here. The archive also has zero-volt baseline
recordings (``_BL`` / ``_Baseline``), ``FreqOnly_*`` and the lab's ``Reduced_*`` outputs,
plus per-recording ``.jpg`` plots.
"""

from __future__ import annotations

import contextlib
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "data/acoustic/US_Data.zip"
LOCAL = ROOT / "data/acoustic"
_TARGET = re.compile(r"/(\d{6})_FUS_AAV_Controller/Mouse_Cntr_(\d+)_(Target[^/]*)\.mat$")


@dataclass(frozen=True)
class Member:
    day: str          # yyyymmdd
    mouse: int
    recording: str    # "Target4", "Target2_Repeat", "Target5Repeat"
    name: str         # path inside the zip
    size: int

    @property
    def filename(self) -> str:
        return f"Mouse_Cntr_{self.mouse:02d}_{self.recording}.mat"


def _day(mmddyy: str) -> str:
    return f"20{mmddyy[4:6]}{mmddyy[0:2]}{mmddyy[2:4]}"


def recordings(archive: Path = ARCHIVE) -> list[Member]:
    """Every sonication recording in the archive, sorted by mouse then name."""
    out = []
    with zipfile.ZipFile(archive) as z:
        for i in z.infolist():
            m = _TARGET.search(i.filename)
            if m:
                out.append(Member(_day(m[1]), int(m[2]), m[3], i.filename, i.file_size))
    return sorted(out, key=lambda r: (r.mouse, r.recording))


@contextlib.contextmanager
def local_path(member: Member, archive: Path = ARCHIVE) -> Iterator[Path]:
    """A readable path for ``member``: the unpacked copy if present, else a temporary one."""
    local = LOCAL / member.day / member.filename
    if local.exists() and local.stat().st_size == member.size:
        yield local
        return
    with tempfile.TemporaryDirectory(prefix="fus_rec_") as d, zipfile.ZipFile(archive) as z:
        path = Path(d) / member.filename
        with z.open(member.name) as src, open(path, "wb") as dst:
            shutil.copyfileobj(src, dst, 16 << 20)
        yield path


def find(mouse: int, recording: str, archive: Path = ARCHIVE) -> Member:
    """The archive member for ``Mouse_Cntr_<mouse>_<recording>.mat``."""
    for m in recordings(archive):
        if m.mouse == mouse and m.recording == recording:
            return m
    raise FileNotFoundError(f"Mouse_Cntr_{mouse:02d}_{recording}.mat is not in {archive.name}")


def from_filename(name: str, archive: Path = ARCHIVE) -> Member:
    """Look a recording up by file name, e.g. ``Mouse_Cntr_01_Target1.mat``."""
    m = re.fullmatch(r"Mouse_Cntr_(\d+)_(Target[^/]*?)(\.mat)?", Path(name).name)
    if not m:
        raise FileNotFoundError(name)
    return find(int(m[1]), m[2], archive)


@contextlib.contextmanager
def paths(mouse: int, recs: list[str], archive: Path = ARCHIVE) -> Iterator[list[Path]]:
    """Readable paths for several recordings of one mouse at once (e.g. a repeated target)."""
    with contextlib.ExitStack() as stack:
        yield [stack.enter_context(local_path(find(mouse, r, archive), archive)) for r in recs]
