"""Turning a profile into the five numbers that get reported.

CaverDock writes an energy for every disc of the tunnel. What a reader wants is where the barrier
is, how much it costs to get in, and how much better the site is than standing outside. Those come
from three points of the profile:

    E_surface   at the mouth of the tunnel
    E_max       the barrier, the highest energy anywhere along it
    E_bound     at the active site

    Ea      = E_max   - E_surface     what entering costs
    dE_BS   = E_bound - E_surface     how much better the destination is

Which end of the profile is the mouth is not something to assume. The direction in the file name
says it, but bulk CaverWeb downloads reuse identifiers and the name can be wrong; the tunnel radius
cannot. A tunnel is narrow at its mouth and opens into the cavity, so the wider end is the binding
site. When the radius disagrees with the file name, both are reported and the row is flagged rather
than quietly resolved.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

# Below this a "tunnel" is the mouth of the pocket rather than a path through the protein. CAVER
# still scores it well -- there is no distance to cross, so nothing obstructs -- which makes it the
# most common way to read a transport result backwards.
SHORT_TUNNEL_A = 2.0


@dataclass
class Metrics:
    energy_surface: Optional[float] = None
    energy_bound: Optional[float] = None
    energy_max: Optional[float] = None
    activation: Optional[float] = None          # Ea
    delta_bs: Optional[float] = None            # dE_BS
    orientation: str = ""                       # what the radii say
    n_discs: int = 0
    span: Optional[float] = None                # how far the profile runs, in angstrom
    flags: tuple = ()


def _energies(profile, bound: str):
    """LB energies, always ordered surface first."""
    values = [p.energy_lb for p in profile]
    return values if bound == "last" else list(reversed(values))


def orientation_from_radius(profile) -> str:
    """Which end of the profile is the binding site, judged by the tunnel's own geometry."""
    if len(profile) < 2:
        return "unknown"
    first, last = profile[0].radius, profile[-1].radius
    if abs(first - last) < 1e-9:
        return "unknown"
    return "last" if last > first else "first"


def evaluate(job, tunnel=None) -> Metrics:
    """The reported numbers for one job, with what could mislead them attached."""
    profile = job.profile
    if not profile:
        return Metrics(flags=("failed",))

    bound = orientation_from_radius(profile)
    flags = []
    if bound == "unknown":
        # A tunnel of constant radius says nothing; fall back to what the file name claims.
        bound = "last" if job.direction == "in" else "first"
        flags.append("orientation_from_name")
    elif job.direction in ("in", "out"):
        expected = "last" if job.direction == "in" else "first"
        if bound != expected:
            flags.append("direction_mismatch")

    values = _energies(profile, bound)
    e_surface, e_bound = values[0], values[-1]
    e_max = max(values)

    if e_surface > 0:
        # Already clashing at the mouth. dE_BS then looks excellent for the wrong reason: it is
        # large because a positive number was subtracted, not because the site is favourable.
        flags.append("positive_surface")
    if not job.has_ub:
        flags.append("lower_bound_only")
    if tunnel is not None and tunnel.length < SHORT_TUNNEL_A:
        flags.append("short_tunnel")

    return Metrics(energy_surface=e_surface, energy_bound=e_bound, energy_max=e_max,
                   activation=e_max - e_surface, delta_bs=e_bound - e_surface,
                   orientation=bound, n_discs=len(profile),
                   span=abs(profile[-1].distance - profile[0].distance),
                   flags=tuple(flags))


def coverage(jobs) -> dict:
    """Which ligand x tunnel x direction combinations are missing, and which arrived twice.

    A CaverWeb combination that fails leaves no log behind, so the only trace of it is the gap.
    Duplicates matter for the opposite reason: two archives claiming the same combination means at
    least one identifier was reused, and the numbers cannot both be right.
    """
    seen, receptors = {}, {}
    for job in jobs:
        if job.tunnel is None or job.direction is None:
            continue
        key = (job.receptor, job.ligand, job.tunnel, job.direction)
        seen.setdefault(key, []).append(job.source)
        r = receptors.setdefault(job.receptor, {"ligands": set(), "tunnels": set()})
        r["ligands"].add(job.ligand)
        r["tunnels"].add(job.tunnel)

    missing, duplicated = [], []
    for receptor, info in receptors.items():
        for ligand in sorted(info["ligands"]):
            for tunnel in sorted(info["tunnels"]):
                for direction in ("in", "out"):
                    key = (receptor, ligand, tunnel, direction)
                    if key not in seen:
                        missing.append(key)
                    elif len(seen[key]) > 1:
                        duplicated.append((key, seen[key]))
    return {"missing": missing, "duplicated": duplicated,
            "expected": sum(len(i["ligands"]) * len(i["tunnels"]) * 2 for i in receptors.values()),
            "present": len(seen)}
