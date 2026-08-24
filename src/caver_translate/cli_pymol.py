"""Write PyMOL scripts for CaverDock trajectories, with the poses chosen by the energy profile.

    caver-pymol CaverWEB/8HTB/met3in4ywxjawqzf_results.zip --object MethylEsterT3In

More than one archive at a time, with the shell choosing the set. The archive names carry the
compound, the tunnel and the direction, so a glob is the whole selection language and there is
nothing new to learn:

    caver-pymol CaverWEB/8HTB/*3in*_results.zip --session pymol_8HTB.pml -o poses/
    caver-pymol CaverWEB/8HTB/benzo*_results.zip --session pymol_8HTB.pml -o poses/
    caver-pymol CaverWEB/8HTB/ --session pymol_8HTB.pml -o poses/        (the lot)

The object name is what the trajectory is called in the loaded session. CaverWeb names it
``traj_ligand_<hash>`` and a renamed script may call it something readable; the hash in the archive
file name is the same hash, so ``--session`` can look the name up instead of being told it. That is
what makes a whole folder possible: one name per archive, none of them typed.
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

# benzo3inuh6v4fpkxc_results.zip -> benzo, tunnel 3, in, hash uh6v4fpkxc. Every archive CaverWeb
# produces is named this way, and it is the only place the tunnel number is written down.
TAIL = re.compile(r"(?P<compound>[A-Za-z]*)(?P<tunnel>\d+)(?P<dir>in|out)(?P<hash>[A-Za-z0-9]+)_results")


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


def collect(given) -> list:
    """Every archive named, with folders expanded.

    Matched by pattern, not by suffix: a second download of the same job arrives as
    ``..._results (1).zip``, which a ``*_results.zip`` glob skips in silence, losing a whole
    calculation from the batch.
    """
    out = []
    for path in map(Path, given):
        if path.is_dir():
            out += sorted(path.glob("*_results*.zip"))
        else:
            out.append(path)
    return out


def for_pymol(path: Path) -> str:
    """The same file, written the way the PyMOL that will read it can open it.

    This runs in WSL and PyMOL runs on Windows, so a /mnt/c/... path pasted after an @ is a file
    PyMOL cannot find. The drive letter form names the same file and both sides accept it.
    """
    parts = path.parts
    if len(parts) > 3 and parts[1] == "mnt" and len(parts[2]) == 1 and parts[2].isalpha():
        return parts[2].upper() + ":/" + "/".join(parts[3:])
    return str(path)


def menu(title, options, multi=True):
    """One question: what there is, and the key to type for it.

    The key is shown beside every line and is what gets typed. Nothing is chosen by a hidden
    position in a list whose order changes with the data.
    """
    print()
    print(title)
    for key, label in options:
        print(f"  {key:>3}   {label}")
    keys = [k for k, _ in options]
    while True:
        raw = input("  > (Enter for all) " if multi else "  > ").strip().lower()
        if not raw and multi:
            return keys
        if raw in keys:
            return [raw]
        print("  Type one of: " + ", ".join(keys))


def folders_with_archives(root: Path) -> list:
    """Where the results are: here, or one folder down.

    One folder down because the archives sit per protein -- CaverWEB/8HTB, CaverWEB/3SQY -- and
    being able to start from either level saves the one cd that is easy to get wrong.
    """
    here = [root] if any(root.glob("*_results*.zip")) else []
    return here + sorted(d for d in root.iterdir()
                         if d.is_dir() and any(d.glob("*_results*.zip")))


def pick(title, paths, label):
    """One of a list of files, or the only one there is."""
    if len(paths) == 1:
        return paths[0]
    keyed = {str(i): p for i, p in enumerate(paths, 1)}
    chosen = menu(title, [(k, label(p)) for k, p in keyed.items()], multi=False)
    return keyed[chosen[0]]


def interactive(args) -> bool:
    """Ask which trajectories to draw, and fill in args as if they had been typed.

    Everything asked about is read off the disk: the archive names carry the compound, the tunnel
    and the direction, so there is nothing to configure and nothing that can go stale.
    """
    folders = folders_with_archives(Path.cwd())
    if not folders:
        print("No CaverDock archives (*_results*.zip) here, or one folder down.", file=sys.stderr)
        return False
    folder = pick("Results folders:", folders, lambda p: p.name)

    # A session script is one that loads trajectories. Chosen by that and not by name, because
    # the scripts this program writes are .pml too and sit in a folder of their own right here:
    # offered as sessions they would renumber the menu the moment the first one existed.
    sessions = [f for f in sorted(folder.glob("*/*.pml")) + sorted(folder.glob("*.pml"))
                if LOAD_LINE.search(f.read_text(errors="ignore"))]
    if not sessions:
        print(f"No PyMOL session script (.pml) under {folder}. The trajectory names are looked",
              file=sys.stderr)
        print("up in it, so one is needed.", file=sys.stderr)
        return False
    session = pick("PyMOL session to take the object names from:", sessions,
                   lambda p: f"{p.parent.name}/{p.name}")

    jobs = [(TAIL.search(a.name), a) for a in sorted(folder.glob("*_results*.zip"))]
    unreadable = [a.name for m, a in jobs if not m]
    jobs = [(m, a) for m, a in jobs if m]
    if not jobs:
        print(f"No archive in {folder} is named the way CaverWeb names them.", file=sys.stderr)
        return False
    if unreadable:
        print(f"Skipping {len(unreadable)} archive(s) not named like CaverWeb's: "
              + ", ".join(unreadable))

    def count(field, value):
        return sum(1 for m, _ in jobs if m[field] == value)

    tunnels = sorted({m["tunnel"] for m, _ in jobs}, key=int)
    names = sorted({m["compound"] for m, _ in jobs})
    keyed = {str(i): c for i, c in enumerate(names, 1)}

    want_t = menu("Tunnel:",
                  [(t, f"tunnel {t}   ({count('tunnel', t)} trajectories)") for t in tunnels])
    want_c = [keyed[k] for k in
              menu("Compound:",
                   [(k, f"{c}   ({count('compound', c)} trajectories)") for k, c in keyed.items()])]
    want_d = menu("Direction:", [("i", "in    -- towards the active site"),
                                 ("o", "out   -- away from it")])

    chosen = [a for m, a in jobs if m["tunnel"] in want_t and m["compound"] in want_c
              and m["dir"][0] in want_d]
    if not chosen:
        print("Nothing matches that combination.", file=sys.stderr)
        return False

    out = folder / "poses"
    # Made now, so that a run that picked exactly one archive still writes into the folder rather
    # than creating a file called poses.
    out.mkdir(parents=True, exist_ok=True)
    args.results = [str(a) for a in chosen]
    args.session = str(session)
    args.out = str(out)

    print()
    print(f"{len(chosen)} of {len(jobs)} trajectories. The same run without the questions:")
    print(f"  caver-pymol {' '.join(sorted(a.name for a in chosen))} "
          f"--session {session} -o {out}")
    return True


def one(results: Path, args, tail):
    """The script for a single archive, or None with the reason already on stderr."""
    job = parse_job(results)[0]
    if not job.profile:
        print(f"ERROR: {results.name} holds no profile ({job.note}).", file=sys.stderr)
        return None

    obj = args.object
    if not obj:
        obj = object_for(args.session, tail["hash"]) if tail else None
        if not obj:
            print(f"ERROR: {results.name}: no trajectory in the session matches its hash.",
                  file=sys.stderr)
            return None

    tunnel = args.tunnel_object
    if tunnel is None:
        tunnel = f"tun_cl_{tail['tunnel']}" if tail else ""

    bound = orientation_from_radius(job.profile)
    if bound == "unknown":
        bound = "last" if job.direction == "in" else "first"

    return script(obj, job.profile, bound=bound, tunnel_obj=tunnel,
                  receptor_obj=args.receptor_object, extra=args.extra,
                  labels=not args.no_labels)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="caver-pymol", description=__doc__.splitlines()[0])
    ap.add_argument("results", nargs="*",
                    help="CaverDock *_results.zip archives, or folders of them; a shell glob "
                         "picks the tunnel and the compounds. Give none and it asks")
    ap.add_argument("--object", help="the trajectory's name in the loaded PyMOL session "
                                     "(one archive only; otherwise use --session)")
    ap.add_argument("--session", help="a pymol .pml to look the names up in, by hash")
    ap.add_argument("--tunnel-object", default=None,
                    help="the tunnel to show, as a mesh; the others are switched off. "
                         "Default: the tun_cl_N named in the archive; pass '' for no tunnel")
    ap.add_argument("--receptor-object", default="structure",
                    help="the protein to fade to 80%% transparent and strip of waters "
                         "(default: %(default)s; empty to leave it alone)")
    ap.add_argument("--extra", type=int, default=0, metavar="N",
                    help="context poses between the three that matter (default: none)")
    ap.add_argument("--no-labels", action="store_true", help="no text beside each pose")
    ap.add_argument("-o", "--out", help="a file, or a folder for one script per archive")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    if not args.results:
        try:
            if not interactive(args):
                return 1
        except (EOFError, KeyboardInterrupt):
            print("Nothing written.", file=sys.stderr)
            return 1

    archives = collect(args.results)
    if not archives:
        print("ERROR: no results archives found.", file=sys.stderr)
        return 1
    missing = [a for a in archives if not a.is_file()]
    if missing:
        print(f"ERROR: not a file: {missing[0]}", file=sys.stderr)
        return 1
    if not args.object and not args.session:
        print("ERROR: give --object, or --session to look the names up.", file=sys.stderr)
        return 1
    if args.object and len(archives) > 1:
        print("ERROR: --object names one trajectory. For more than one archive use --session.",
              file=sys.stderr)
        return 1

    out = Path(args.out).resolve() if args.out else None
    if len(archives) > 1 and out is None:
        print("ERROR: more than one archive needs -o, a folder to write them into.",
              file=sys.stderr)
        return 1
    folder = out is not None and (len(archives) > 1 or out.is_dir())
    if folder:
        out.mkdir(parents=True, exist_ok=True)

    failed = 0
    for results in archives:
        text = one(results, args, TAIL.search(results.name))
        if text is None:
            failed += 1
            continue
        if out is None:
            print(text, end="")
            continue
        path = out / f"{results.name.split('_results')[0]}.pml" if folder else out
        path.write_text(text, encoding="utf-8")
        # Absolute, because the CaverWeb session script ends with the working directory inside
        # data/ and never comes back: a relative @path is looked for there and is not found.
        print(f"{results.name}  ->  @{for_pymol(path)}")

    if failed:
        print(f"{failed} of {len(archives)} archives produced nothing.", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
