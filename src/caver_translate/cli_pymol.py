"""Write a PyMOL script for one trajectory, with the poses chosen by the energy profile.

    caver-pymol CaverWEB/8HTB/met3in4ywxjawqzf_results.zip --object MethylEsterT3In

The object name is what the trajectory is called in the loaded session. CaverWeb names it
``traj_ligand_<hash>`` and a renamed script may call it something readable; the hash in the archive
file name is the same hash, so ``--session`` can look the name up instead of being told it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from . import __version__
from .metrics import orientation_from_radius
from .parse import parse_job
from .figures import script

LOAD_LINE = re.compile(r'load\s+"?trajectory/ligand_(?P<hash>[A-Za-z0-9]+)\.pdbqt"?\s*,\s*"?(?P<obj>[^"\s]+)"?')


def object_for(session_pml, hash_: str):
    """The name this trajectory has in a session script, found by its hash.

    CaverWeb writes ``load "trajectory/ligand_<hash>.pdbqt", "traj_ligand_<hash>"``. Renaming the
    objects to something readable rewrites the second half and leaves the first, so the hash is the
    one thing that survives and it is the same hash the results archive is named after.
    """
    for line in Path(session_pml).read_text(errors="ignore").splitlines():
        m = LOAD_LINE.search(line)
        if m and m.group("hash") == hash_:
            return m.group("obj")
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="caver-pymol", description=__doc__.splitlines()[0])
    ap.add_argument("results", help="a CaverDock *_results.zip")
    ap.add_argument("--object", help="the trajectory's name in the loaded PyMOL session")
    ap.add_argument("--session", help="a pymol .pml to look the name up in, by hash")
    ap.add_argument("--tunnel-object", default="",
                    help="the tunnel to show, as a mesh; the others are switched off, e.g. tun_cl_3")
    ap.add_argument("--receptor-object", default="structure",
                    help="the protein to fade to 80%% transparent and strip of waters "
                         "(default: %(default)s; empty to leave it alone)")
    ap.add_argument("--extra", type=int, default=0, metavar="N",
                    help="context poses between the three that matter (default: none)")
    ap.add_argument("--no-labels", action="store_true", help="no text beside each pose")
    ap.add_argument("-o", "--out", help="write here instead of to the screen")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    results = Path(args.results)
    if not results.is_file():
        print(f"ERROR: {results} is not a file.", file=sys.stderr)
        return 1

    jobs = parse_job(results)
    job = jobs[0]
    if not job.profile:
        print(f"ERROR: {results.name} holds no profile ({job.note}).", file=sys.stderr)
        return 1

    obj = args.object
    if not obj and args.session:
        m = re.search(r"(?:in|out)([A-Za-z0-9]+)_results", results.name)
        if m:
            obj = object_for(args.session, m.group(1))
        if not obj:
            print(f"ERROR: no trajectory in {args.session} matches this archive's hash.",
                  file=sys.stderr)
            return 1
    if not obj:
        print("ERROR: give --object, or --session to look it up.", file=sys.stderr)
        return 1

    bound = orientation_from_radius(job.profile)
    if bound == "unknown":
        bound = "last" if job.direction == "in" else "first"

    text = script(obj, job.profile, bound=bound, tunnel_obj=args.tunnel_object,
                  receptor_obj=args.receptor_object, extra=args.extra,
                  labels=not args.no_labels)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Written to {args.out}. In PyMOL, with the session already loaded:  @{args.out}")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
