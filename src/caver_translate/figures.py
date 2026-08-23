"""Pick the poses worth drawing, and write the PyMOL script that draws them.

Named figures rather than pymol on purpose: a module called pymol inside this package
shadows the real PyMOL for anything in it that imports the name.

A CaverDock trajectory is one pose per disc of the tunnel -- sixty-eight of them for a thirteen
angstrom route. Taking five at even spacing is easy and says nothing: the pose that matters is the
one at the top of the energy profile, and even spacing lands on it only by luck.

The states in the PyMOL object are the points of the profile, in order and one for one (verified on
real trajectories: 68 discs, 68 MODEL records). So the profile can choose the states, and each
choice can carry its reason.

What comes out is a plain .pml with the state numbers already resolved and a comment beside each
saying why it was picked. Editing it needs no Python and no help: the numbers are right there.

Seven-bit ASCII throughout. PyMOL read a middle dot in a label as latin-1 and printed mojibake.
"""
from __future__ import annotations

from pathlib import Path

# Blue to orange, ordered so that progress along the tunnel reads left to right, and none of them
# is red -- the tunnel mesh is usually red, and a red ligand on it disappears.
PALETTE = ["marine", "cyan", "green", "yellow", "orange", "magenta", "salmon", "purple"]

# Two poses within this of each other are the same picture twice. Vina reports to one decimal,
# so anything under a few tenths is not a difference a figure can show.
MEANINGFUL_KCAL = 0.3


def _fmt(value, digits=1) -> str:
    return "?" if value is None else f"{value:.{digits}f}"


def choose_states(profile, bound: str = "last", extra: int = 0) -> list:
    """The states worth showing, as (state, tag, label, reason).

    Three always: where it starts, the barrier, where it ends. The barrier is the point of the
    exercise -- it is what a reader is being shown -- and it is the one an evenly spaced sample
    misses. ``extra`` adds evenly spaced context poses between them for a fuller figure.

    ``state`` is one-based, the way PyMOL counts. ``label`` is what goes on screen: distance from
    the active site and the energy there, because those are the two numbers a figure is read for.
    """
    if not profile:
        return []
    n = len(profile)
    energies = [p.energy_lb for p in profile]

    surface_i = 0 if bound == "last" else n - 1
    bound_i = n - 1 if bound == "last" else 0
    barrier_i = max(range(n), key=lambda i: energies[i])
    lowest_i = min(range(n), key=lambda i: energies[i])

    picked = {
        surface_i: ("start", "entrance", "at the mouth of the tunnel"),
        barrier_i: ("barrier", "HARDEST STEP", "highest energy on the profile: the obstacle"),
        bound_i: ("end", "active site", "at the binding site"),
    }
    # The deepest point is usually the bound end itself. Adding it then draws a second pose on
    # top of the first for no gain, so it earns a place only by being meaningfully deeper.
    if energies[lowest_i] < energies[bound_i] - MEANINGFUL_KCAL:
        picked.setdefault(lowest_i, ("lowest", "deepest point",
                                     "lowest energy on the profile, deeper than the site"))

    for k in range(extra):
        i = round((k + 1) * (n - 1) / (extra + 1))
        picked.setdefault(i, ("step", "", "context, evenly spaced"))

    # Distance is measured from the active site, not from wherever the file happens to start: that
    # is the axis the figure is read along, and for a run stored backwards the two are opposite.
    site = profile[bound_i].distance

    out = []
    for i in sorted(picked):
        tag, title, reason = picked[i]
        point = profile[i]
        to_site = abs(site - point.distance)
        text = f"{_fmt(to_site)} A from site   {_fmt(point.energy_lb, 1)} kcal/mol"
        if title:
            text = f"{title}   {text}"
        out.append((i + 1, tag, text,
                    f"{reason} | disc {point.disc}, {_fmt(to_site)} A from the site, "
                    f"{_fmt(point.energy_lb, 2)} kcal/mol"))
    return out


def script(obj: str, profile, bound: str = "last", tunnel_obj: str = "",
           receptor_obj: str = "structure", extra: int = 0, labels: bool = True,
           prefix: str = "snap") -> str:
    """A ready-to-run .pml for one trajectory."""
    states = choose_states(profile, bound=bound, extra=extra)
    if not states:
        return "# no profile: nothing to draw\n"

    needed = [n for n in (obj, tunnel_obj) if n]
    L = [
        "# ---------------------------------------------------------------------------",
        f"# Poses of {obj}, chosen from its energy profile.",
        "#",
        "# Each state below is a point of the profile, and the comment says why it was picked.",
        "# To use fewer, delete lines. To use others, change the number: state N is profile point",
        "# N, counting from the mouth of the tunnel.",
        "#",
        "# Colours, in order of appearance: " + ", ".join(PALETTE[:len(states)]),
        "# ---------------------------------------------------------------------------",
        "",
        "# Does the session this needs actually exist? Without it every line below fails on its",
        "# own and the screen fills with errors that each describe a symptom, never the cause.",
        "python",
        "from pymol import cmd",
        f"_needed = {needed!r}",
        "_missing = [n for n in _needed if n not in cmd.get_object_list()]",
        "if _missing:",
        "    print('')",
        "    print('  This script draws poses from a session that is not loaded.')",
        "    print('  Missing: ' + ', '.join(_missing))",
        "    print('')",
        "    print('  Load the CaverWeb session first, from the folder that contains data/:')",
        "    print('     cd /path/to/pymol_<jobid>')",
        "    print('     @pymol.pml            (or your renamed copy)')",
        "    print('  then run this script again.')",
        "    print('')",
        "python end",
        "",
        f"delete {prefix}_*",
        "delete pose_label_*",
        "set all_states, on",
        "",
    ]

    if receptor_obj:
        L += [
            "# The receptor, faint: the figure is of what happens inside it, and a solid protein",
            "# hides the tunnel and every pose in it. Waters go because they are not part of the",
            "# story and each one is another dot in front of the route.",
            f"show cartoon, {receptor_obj}",
            f"hide lines, {receptor_obj}",
            f"hide spheres, {receptor_obj}",
            f"hide nonbonded, {receptor_obj}",
            f"remove {receptor_obj} and solvent",
            f"set cartoon_transparency, 0.8, {receptor_obj}",
            f"color grey80, {receptor_obj}",
            "",
        ]

    if tunnel_obj:
        L += [
            "# Only the tunnel this trajectory goes through. The others are still loaded and can be",
            "# switched back on with: enable tun_cl_*",
            "hide everything, tun_cl_*",
            "disable tun_cl_*",
            f"enable {tunnel_obj}",
            f"show mesh, {tunnel_obj}",
            "set mesh_width, 0.3",
            f"color grey70, {tunnel_obj}",
            "",
        ]

    for i, (state, tag, text, reason) in enumerate(states):
        colour = PALETTE[i % len(PALETTE)]
        name = f"{prefix}_{i + 1}_{tag}"
        L += [
            f"# {reason}",
            f"create {name}, {obj}, {state}, 1",
            f"show sticks, {name}",
            f"color {colour}, {name} and elem C",
        ]
        if labels:
            L.append(f'pseudoatom pose_label_{i + 1}, selection=({name}), label="{text}"')
        L.append("")

    L += [
        "# Text rather than markers: a sphere at each end is one more object in front of the route,",
        "# and the labels already say which end is which.",
        "set label_size, 15",
        "set label_color, white",
        "set label_outline_color, black",
        "set label_position, (0, 1.8, 0)",
        "",
        "# orient on the route, then:  ray 1600, 1200   and   png figure.png, dpi=300",
        f"orient {prefix}_*",
    ]
    return "\n".join(L) + "\n"


def write_script(path, obj: str, profile, **kwargs) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(script(obj, profile, **kwargs), encoding="utf-8")
    return path
