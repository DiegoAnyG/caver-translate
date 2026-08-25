"""Point it at CaverDock output and get back something readable.

    caver-translate <folder> -o out/

A CaverWeb download -- eighty archives named after hashes -- or a run made on this machine, which
CaverWeb's summary file never reaches. Either becomes two CSVs and one page. Nothing is modified in
the folder it reads.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .metrics import coverage
from .parse import scan
from .report import COLUMNS, TUNNEL_COLUMNS, rows, tunnel_rows, write_csv, write_html


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="caver-translate", description=__doc__.splitlines()[0])
    ap.add_argument("folder", help="a CaverWeb download (one sub-folder per receptor, or a single "
                                   "receptor's folder), or a folder of CaverDock output produced "
                                   "here")
    ap.add_argument("-o", "--out", default="caver_report", help="where to write (default: %(default)s)")
    ap.add_argument("--no-html", action="store_true", help="tables only")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    folder = Path(args.folder)
    if not folder.is_dir():
        print(f"ERROR: {folder} is not a folder.", file=sys.stderr)
        return 1

    tunnels, jobs = scan(folder)
    if not jobs:
        print(f"Nothing to read in {folder}. Expected either a CaverWeb download -- one sub-folder "
              "per receptor, each holding *_results.zip -- or CaverDock output produced here, "
              "which is a folder with a *-lb.pdbqt trajectory or a profile .dat in it.",
              file=sys.stderr)
        return 1

    out = Path(args.out)
    records = rows(tunnels, jobs)
    write_csv(out / "transport.csv", records, COLUMNS)
    if tunnels:
        write_csv(out / "tunnels.csv", tunnel_rows(tunnels), TUNNEL_COLUMNS)
    if not args.no_html:
        write_html(out / "report.html", tunnels, jobs)

    cov = coverage(jobs)
    failed = sum(1 for r in records if "failed" in r["flags"])
    print(f"{len(records)} calculations, {len(tunnels)} tunnels -> {out}/")
    print(f"  combinations present : {cov['present']} of {cov['expected']}")
    if cov["missing"]:
        print(f"  missing              : {len(cov['missing'])} (listed in report.html)")
    if cov["duplicated"]:
        print(f"  claimed twice        : {len(cov['duplicated'])}")
    if failed:
        print(f"  archives with no profile: {failed}")
    flagged = sorted({f for r in records for f in r["flags"].split() if f != "failed"})
    if flagged:
        print("  worth reading before quoting a number: " + ", ".join(flagged))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
