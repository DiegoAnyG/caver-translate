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
    """The states worth showing, as (state, tag, reason).

    Three always: where it starts, the barrier, where it ends. The barrier is the point of the
    exercise -- it is what a reader is being shown -- and it is the one an evenly spaced sample
    misses. ``extra`` adds evenly spaced context poses between them for a fuller figure.

    ``state`` is one-based, the way PyMOL counts.
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
        surface_i: ("start", "at the mouth of the tunnel"),
        barrier_i: ("barrier", "highest energy on the profile: the obstacle"),
        bound_i: ("end", "at the binding site"),
    }
    # The deepest point is usually the bound end itself. Adding it then draws a second pose on
    # top of the first for no gain, so it earns a place only by being meaningfully deeper.
    if energies[lowest_i] < energies[bound_i] - MEANINGFUL_KCAL:
        picked.setdefault(lowest_i, ("lowest", "lowest energy on the profile, deeper than the site"))

    for k in range(extra):
        i = round((k + 1) * (n - 1) / (extra + 1))
        picked.setdefault(i, ("step", "context, evenly spaced"))

    out = []
    for i in sorted(picked):
        tag, reason = picked[i]
        point = profile[i]
        out.append((i + 1, tag,
                    f"{reason} | disc {point.disc}, {_fmt(point.distance)} A, "
                    f"{_fmt(point.energy_lb, 2)} kcal/mol"))
    return out


def script(obj: str, profile, bound: str = "last", tunnel_obj: str = "", extra: int = 0,
           labels: bool = True, prefix: str = "snap") -> str:
    """A ready-to-run .pml for one trajectory."""
    states = choose_states(profile, bound=bound, extra=extra)
    if not states:
        return "# no profile: nothing to draw\n"

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
        f"delete {prefix}_*",
        "delete pose_label_*",
        "delete route_start",
        "delete route_end",
        "set all_states, on",
        "",
    ]

    if tunnel_obj:
        L += [
            f"# The tunnel, as a mesh rather than a surface: a solid tunnel hides the ligand inside it.",
            f"hide surface, {tunnel_obj}",
            f"show mesh, {tunnel_obj}",
            "set mesh_width, 0.3",
            f"color grey70, {tunnel_obj}",
            "",
        ]

    for i, (state, tag, reason) in enumerate(states):
        colour = PALETTE[i % len(PALETTE)]
        name = f"{prefix}_{i + 1}_{tag}"
        L += [
            f"# {reason}",
            f"create {name}, {obj}, {state}, 1",
            f"show sticks, {name}",
            f"color {colour}, {name} and elem C",
        ]
        if labels:
            text = f"{tag} · {reason.split('| ')[-1]}"
            L.append(f'pseudoatom pose_label_{i + 1}, selection=({name}), label="{text}"')
        L.append("")

    first, last = states[0], states[-1]
    L += [
        "# Where the route begins and ends, as spheres, so the direction is readable in a still.",
        f"pseudoatom route_start, selection=({prefix}_1_{first[1]}), label=\"START\"",
        f"pseudoatom route_end, selection=({prefix}_{len(states)}_{last[1]}), label=\"END\"",
        "show spheres, route_start route_end",
        "set sphere_scale, 0.45, route_start",
        "set sphere_scale, 0.45, route_end",
        "color white, route_start",
        "color black, route_end",
        "",
        "set label_size, 16",
        "set label_color, white",
        "set label_outline_color, black",
        "set label_position, (0, 1.6, 0)",
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
