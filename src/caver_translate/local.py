"""Reading CaverDock output produced on this machine.

A CaverWeb download hands back ``results.json``, a summary CaverWeb writes itself. CaverDock does
not produce that file, so a job run locally has to be read from what CaverDock actually leaves
behind, under different names:

    <name>-lb.pdbqt     the lower-bound trajectory, one MODEL per disc
    <name>-ub.pdbqt     the upper-bound trajectory
    <name>-lb.dat       the resolved lower-bound energy profile, if it was written
    <name>-ub.dat       the resolved upper-bound energy profile
    results/profile.dat the same thing under cd-screening, one bound per screening
    *.dsd               the discretised tunnel: centre, normal and radius per disc

``cd-analysis`` and ``cd-screening`` write the ``.dat`` files; a bare ``caverdock`` run writes only
the trajectories, and ``cd-energyprofile -d tunnel.dsd -t <name>-lb.pdbqt`` turns one into the
other. All of them carry the same six columns, which is CaverDock's own reader in
``pycaverdock/energy_profile.py``:

    distance  disc  minE  maxE  radius  lbE

**The energies and radii are taken from the trajectory, not from those columns.** Both agree on a
``cd-analysis`` result -- all 68 discs of the reference job, exactly -- but they do not always, and
where they differ the trajectory is right:

- In a ``cd-screening`` ``profile.dat`` the ``radius`` and ``lbE`` columns hold one value repeated
  down every row, and ``minE`` comes out as ``min(trajectory, lbE)``. That clips the profile at the
  free-docking energy exactly where it rises -- at the barrier. On the run measured here it moved
  ``E_max`` from -2.0 to -3.2 and ``Ea`` from 2.6 to 1.4, and cd-screening's own ``results.csv``,
  built from the same columns, reported ``E_bound = E_max = E_surface`` and ``dE_BS = 0``.
- The distance column is sound in every file seen, so it is what the ``.dat`` is read for. Without
  one the distances are measured off the discretised tunnel instead.

The upper bound is the one thing the trajectory cannot give: it is the search, not the result. It
revisits discs, and neither the first nor the highest energy recorded for one matches the resolved
profile (59 of 68 discs disagree). Upper-bound numbers therefore come from a ``.dat`` and from
nowhere else; without one the job is reported as lower-bound only, which is true.
"""
from __future__ import annotations

import re
from pathlib import Path

from .parse import JOB_NAME, Job, Point

# The first of these in a MODEL is that model's disc: <disc> <energy> <radius>. A trajectory may
# carry a second one repeating disc 0 as a reference, which is why the first is taken and the rest
# of the MODEL ignored.
TUNNEL_REMARK = re.compile(r"^REMARK CAVERDOCK TUNNEL:\s+(\d+)\s+(\S+)\s+(\S+)")

# cd-screening names each experiment folder r<receptor>-l<ligand>-t<tunnel>-d<direction>-<bound>,
# which is the only place a local batch records what it calculated.
EXPERIMENT_NAME = re.compile(r"^r(?P<receptor>.*?)-l(?P<ligand>.*?)-t(?P<tunnel>.*?)"
                             r"-d(?P<direction>in|out)-(?P<bound>lowerbound|upperbound)$")


def _numbers(path) -> list:
    """Rows of floats, skipping comments and anything that is not one."""
    rows = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.lstrip().startswith("#"):
            continue
        try:
            rows.append([float(f) for f in line.split()])
        except ValueError:
            continue
    return rows


def travelled(discs) -> list:
    """How far along the tunnel each disc sits, from the rows of a .dsd.

    Not the centre-to-centre distance: CaverDock measures each step along the disc normal, which is
    shorter wherever the tunnel bends. Against a CaverWeb profile of the same tunnel the projection
    lands within 0.02 % end to end, where the plain sum of centre distances runs 3.5 % long.
    """
    out, total = [0.0] * bool(discs), 0.0
    for a, b in zip(discs, discs[1:]):
        total += abs(sum((b[i] - a[i]) * a[3 + i] for i in range(3)))
        out.append(total)
    return out


def disc_distances(dsd_path) -> list:
    """How far along the tunnel each disc of this .dsd sits."""
    return travelled([r for r in _numbers(dsd_path) if len(r) >= 7])


def trajectory(pdbqt_path) -> list:
    """(disc, energy, radius) per MODEL of a CaverDock trajectory."""
    out, done = [], True
    for line in Path(pdbqt_path).read_text(errors="ignore").splitlines():
        if line.startswith("MODEL"):
            done = False
        elif not done:
            m = TUNNEL_REMARK.match(line)
            if m:
                out.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
                done = True
    return out


def discretised(folder: Path, discs: int) -> list:
    """The tunnel this profile was docked into, wherever in the folder it was written.

    A run can leave more than one: cd-screening keeps the discretised tunnel and then the extended
    one it actually used. The one with a disc for every point of the profile is that one.
    """
    for path in sorted(Path(folder).rglob("*.dsd")):
        rows = [r for r in _numbers(path) if len(r) >= 7]
        if len(rows) == discs:
            return rows
    return []


def _identify(folder: Path) -> tuple:
    """Receptor, ligand, tunnel and direction, from whatever the folder name says.

    Nothing inside a CaverDock output names the calculation, so the folder is the only witness.
    A cd-screening experiment folder names all four; anything else contributes what it can, and the
    rest stays None so the coverage check reports it as unknown rather than inventing it.
    """
    m = EXPERIMENT_NAME.match(folder.name)
    if m:
        digits = re.search(r"\d+", m.group("tunnel"))
        return (m.group("receptor"), m.group("ligand"),
                int(digits.group()) if digits else None, m.group("direction"))
    m = JOB_NAME.match(folder.name)
    if m:
        return folder.parent.name, m.group("ligand"), int(m.group("tunnel")), m.group("direction")
    return folder.parent.name, folder.name, None, None


def _nearest(folder: Path, pattern: str):
    """The shallowest match, so a folder's own files win over a run kept inside it."""
    found = sorted(folder.rglob(pattern), key=lambda p: (len(p.parts), p.name))
    return found[0] if found else None


def _sources(folder: Path) -> tuple:
    """The lower-bound trajectory and the two profile files, wherever this run wrote them.

    cd-analysis names the profiles after the bound. cd-screening runs one bound per screening and
    writes a single results/profile.dat, so which one it holds is in the experiment folder name.
    """
    lb = _nearest(folder, "*-lb.dat")
    ub = _nearest(folder, "*-ub.dat")
    if lb is None and ub is None:
        screening = folder / "results" / "profile.dat"
        if screening.is_file():
            m = EXPERIMENT_NAME.match(folder.name)
            if m is not None and m.group("bound") == "upperbound":
                ub = screening
            else:
                lb = screening
    return _nearest(folder, "*-lb.pdbqt"), lb, ub


def _distances(path) -> list:
    """The distance column, the one part of a profile file that is sound in every run seen."""
    return [r[0] for r in _numbers(path) if len(r) >= 2]


def _repair_radii(folder: Path, profile: list) -> str:
    """Put the tunnel's own radii back when the profile file repeated one down the column."""
    if len({p.radius for p in profile}) > 1:
        return ""
    discs = discretised(folder, len(profile))
    if not discs:
        return "the radius is the same on every disc and there is no tunnel to check it against"
    for point, disc in zip(profile, discs):
        point.radius = round(disc[6], 1)
    return ""


def parse_local_job(folder) -> Job:
    """One local CaverDock output folder as a job."""
    folder = Path(folder)
    receptor, ligand, tunnel, direction = _identify(folder)
    traj, lb_dat, ub_dat = _sources(folder)
    notes = []

    steps = trajectory(traj) if traj is not None else []
    if steps:
        # Positional, not by disc number: a profile file can misnumber its last disc, and does --
        # CaverWeb's results.json calls disc 67 disc 68, and a cd-screening profile.dat skips 78
        # and ends at 79. One row per model is the check that they line up at all.
        along = _distances(lb_dat or ub_dat) if (lb_dat or ub_dat) else []
        if len(along) != len(steps):
            along = travelled(discretised(folder, len(steps)))
            if not along:
                notes.append("no discretised tunnel beside the trajectory, so distances are "
                             "unknown")
        profile = [Point(distance=along[i] if i < len(along) else 0.0, disc=disc, radius=radius,
                         energy_lb=energy)
                   for i, (disc, energy, radius) in enumerate(steps)]
        source = traj.name
    elif lb_dat is not None or ub_dat is not None:
        # No trajectory left on disk. The profile file is all there is, and its lower bound may
        # have been clipped where it rose -- which is where the barrier is.
        path = lb_dat or ub_dat
        energy = 2 if lb_dat is not None else 5           # minE, or the ub file's own lbE column
        profile = [Point(distance=r[0], disc=int(r[1]), radius=r[4], energy_lb=r[energy])
                   for r in _numbers(path) if len(r) > energy]
        source = path.name
        notes.append("read from the profile file, with no trajectory to check it against: a "
                     "barrier clipped at the free-docking energy would not show")
        notes.append(_repair_radii(folder, profile))
    else:
        return Job(receptor, ligand, tunnel, direction, folder.name,
                   note="no CaverDock trajectory or energy profile in this folder")

    if not profile:
        return Job(receptor, ligand, tunnel, direction, folder.name,
                   note=f"{source} held no profile")

    has_ub = ub_dat is not None
    if has_ub:
        bounds = {int(r[1]): (r[2], r[3]) for r in _numbers(ub_dat) if len(r) >= 4}
        for point in profile:
            if point.disc in bounds:
                point.energy_ub_min, point.energy_ub_max = bounds[point.disc]
    else:
        notes.append("no upper-bound profile: run cd-energyprofile on the -ub trajectory to "
                     "write one")
    return Job(receptor=receptor, ligand=ligand, tunnel=tunnel, direction=direction, source=source,
               profile=profile, has_ub=has_ub, note="; ".join(n for n in notes if n))


def local_jobs(folder) -> list:
    """Every local CaverDock output folder under this one.

    Found by the files rather than by the folder layout, so a single job, a cd-screening output
    tree and a folder someone arranged by hand all read the same way. A folder inside one that is
    already a job is part of it -- cd-screening keeps the raw trajectory in an ``intermediate``
    sub-folder -- and is not counted twice.
    """
    root = Path(folder)
    found = {p.parent for p in root.rglob("*-lb.pdbqt")}
    found |= {p.parent for p in root.rglob("*-lb.dat")}
    found |= {p.parent for p in root.rglob("*-ub.dat")}
    found |= {p.parent.parent for p in root.rglob("results/profile.dat")}
    outermost = [d for d in sorted(found) if not any(p in found for p in d.parents)]
    return [parse_local_job(d) for d in outermost]
