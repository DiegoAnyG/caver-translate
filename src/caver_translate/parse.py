"""Reading what CaverWeb hands back.

A CaverWeb download is a folder of zip files whose names are hashes. Each one holds the transport
of one ligand through one tunnel in one direction, and nothing in the file name says which -- so a
receptor with five compounds, three tunnels and two directions arrives as thirty opaque archives
that have to be opened one at a time to find out what they are.

Two things are read here:

- ``<hash>_summary.txt``, the CAVER tunnel table: bottleneck radius, length, curvature, priority.
- ``<hash>_results.zip``, one CaverDock job: ``results.json`` carries the full energy profile,
  point by point, and is the only file needed. The rest of the archive is poses and logs.
"""
from __future__ import annotations

import json
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

# CaverWeb names an individually exported job <ligand><tunnel><direction><hash>. The bulk download
# does not, which is why identity is confirmed against the profile rather than trusted from here.
JOB_NAME = re.compile(r"^(?P<ligand>[A-Za-z-]+?)(?P<tunnel>[0-9]+)(?P<direction>in|out)(?P<hash>[A-Za-z0-9]*)$")

# The CAVER table starts at the line whose first column header is ID; everything above is a legend.
TUNNEL_HEADER = re.compile(r"^\s*ID\s+No\s+No_snaps")

# "..._results.zip", and the copy a browser names "..._results (1).zip".
SUFFIX = re.compile(r"_results(?:\s*\(\d+\))?\.zip$")


@dataclass(frozen=True)
class Tunnel:
    """One tunnel cluster as CAVER reports it."""
    receptor: str
    tunnel: int
    bottleneck_radius: float
    length: float
    curvature: float
    priority: float


@dataclass
class Point:
    """One disc of the discretised tunnel."""
    distance: float
    disc: int
    radius: float
    energy_lb: float
    energy_ub_min: Optional[float] = None
    energy_ub_max: Optional[float] = None


@dataclass
class Job:
    """One CaverDock calculation: this ligand, this tunnel, this direction."""
    receptor: str
    ligand: str
    tunnel: Optional[int]
    direction: Optional[str]
    source: str
    profile: list = field(default_factory=list)
    has_ub: bool = False
    note: str = ""

    @property
    def combo(self) -> str:
        return f"{self.ligand}/t{self.tunnel}/{self.direction}"


def parse_tunnels(summary_path) -> list:
    """The tunnel table from a CAVER summary file.

    Fixed-width columns under a header line, preceded by a legend explaining each one. Rows are
    read by position after the header rather than by splitting the legend, which contains the same
    words and would match.
    """
    path = Path(summary_path)
    receptor = path.parent.name
    out, in_table = [], False
    for line in path.read_text(errors="ignore").splitlines():
        if TUNNEL_HEADER.match(line):
            in_table = True
            continue
        if not in_table:
            continue
        fields = line.split()
        if not fields or not fields[0].isdigit():
            if out:
                break            # a blank line or the citation footer closes the table
            continue
        try:
            out.append(Tunnel(receptor=receptor, tunnel=int(fields[0]),
                              bottleneck_radius=float(fields[3]), length=float(fields[6]),
                              curvature=float(fields[8]), priority=float(fields[10])))
        except (IndexError, ValueError):
            continue
    return out


def _points(profile) -> list:
    return [Point(distance=float(p["distance"]), disc=int(p["disc"]), radius=float(p["radius"]),
                  energy_lb=float(p["energyLb"]),
                  energy_ub_min=p.get("energyUbMin"), energy_ub_max=p.get("energyUbMax"))
            for p in profile]


def parse_job(zip_path) -> list:
    """The jobs inside one results archive.

    Returns a job even when the archive holds no profile: a CaverWeb combination that failed leaves
    no log to inspect, so the only way to report it is to notice the gap where it should have been.
    """
    path = Path(zip_path)
    receptor = path.parent.name
    # A second download of the same archive is "..._results (1).zip": the suffix is not at
    # the end, so trimming a fixed number of characters eats part of the identifier.
    stem = SUFFIX.sub("", path.name)
    m = JOB_NAME.match(stem)
    ligand = m.group("ligand") if m else stem
    tunnel = int(m.group("tunnel")) if m else None
    direction = m.group("direction") if m else None

    def blank(note):
        return [Job(receptor, ligand, tunnel, direction, path.name, note=note)]

    try:
        with zipfile.ZipFile(path) as zf:
            name = next((n for n in zf.namelist() if n.endswith("results.json")), None)
            if name is None:
                return blank("no results.json: the calculation produced nothing")
            data = json.loads(zf.read(name))
    except (zipfile.BadZipFile, json.JSONDecodeError, OSError) as e:
        return blank(f"unreadable: {e}")

    jobs = []
    for entry in data:
        profile = _points(entry.get("profile") or [])
        jobs.append(Job(receptor=receptor,
                        ligand=ligand if len(data) == 1 else f"{ligand}:{entry.get('name')}",
                        tunnel=tunnel, direction=direction, source=path.name,
                        profile=profile, has_ub=bool(entry.get("hasUb")),
                        note="" if profile else "empty profile"))
    return jobs or blank("results.json held no entries")


def result_archives(folder) -> list:
    """The result archives in one folder.

    Matched by pattern rather than by suffix: a second download is named "..._results (1).zip",
    which does not end in "_results.zip" and was silently skipped -- the whole calculation vanished
    from the report, which is the one failure a tool for auditing these downloads must not have.
    The same filter keeps out the other archives CaverWeb ships, the CAVER output and the PyMOL
    session.
    """
    return sorted(p for p in Path(folder).glob("*.zip") if SUFFIX.search(p.name))


def scan(folder) -> tuple:
    """Every tunnel table and every job under a CaverWeb download folder.

    Receptors are the sub-folders: that is how CaverWeb downloads are kept once more than one
    target is involved, and the folder name is the only place the target's name survives.
    """
    root = Path(folder)
    tunnels, jobs = [], []
    for summary in sorted(root.glob("*/*_summary.txt")):
        tunnels.extend(parse_tunnels(summary))
    for receptor_dir in sorted(d for d in root.iterdir() if d.is_dir()):
        for zip_path in result_archives(receptor_dir):
            jobs.extend(parse_job(zip_path))
    if not tunnels and not jobs:                      # a single receptor, given directly
        for summary in sorted(root.glob("*_summary.txt")):
            tunnels.extend(parse_tunnels(summary))
        for zip_path in result_archives(root):
            jobs.extend(parse_job(zip_path))
    return tunnels, jobs
