# caver-translate

Turn a CaverWeb download into a table you can read.

CaverWeb runs one job per **ligand × tunnel × direction** and hands each one back as a zip named
after a hash. Five compounds, three tunnels, two directions is thirty archives whose names say
nothing about what is inside. A modest study is eighty.

This reads the folder and writes what gets reported: what entering costs, where the barrier is,
how good the destination is — with the things that mislead those numbers marked on the rows they
affect.

```bash
pip install caver-translate
caver-translate CaverWEB/ -o report/
```

```
84 calculations, 18 tunnels -> report/
  combinations present : 80 of 90
  missing              : 10 (listed in report.html)
  claimed twice        : 1
  archives with no profile: 3
  worth reading before quoting a number: direction_mismatch, positive_surface, short_tunnel
```

Three files come out: `transport.csv`, `tunnels.csv`, and `report.html` — one page, every profile
drawn, easiest route first. No dependencies: the plots are SVG written by hand, so the report is a
single file that opens anywhere.

## What it reports

| | |
|---|---|
| `E_surface` | energy at the mouth of the tunnel |
| `E_max` | the barrier, the highest energy anywhere along it |
| `E_bound` | energy at the active site |
| `Ea` | `E_max − E_surface` — what entering costs |
| `dE_BS` | `E_bound − E_surface` — how much better the destination is |

Alongside them, the tunnel's own geometry from CAVER: bottleneck radius, length, curvature,
priority. `Ea` without the length beside it is how a tunnel that is barely a tunnel wins.

## What it flags

Each of these was met in real data before it was written down.

- **`short_tunnel`** — under 2 Å. CAVER gives it the best priority because there is no distance in
  which to be obstructed. It is the mouth of the pocket, not a route through the protein.
- **`positive_surface`** — the ligand already clashes at the mouth, so `dE_BS` looks excellent
  because a positive number was subtracted, not because the site binds well.
- **`direction_mismatch`** — the tunnel widens the opposite way to what the file name claims. A
  tunnel is narrow at its mouth and opens into the cavity, so the radius decides and the name is
  reported as suspect. This found a genuinely mislabelled archive on its first run.
- **`lower_bound_only`** — no upper-bound trajectory. The lower bound can pass through
  discontinuities and understate a barrier.
- **`failed`** — the archive holds no profile. CaverWeb writes no log for a combination that
  fails, so the empty archive and the gap in the table are the only evidence it was attempted.

Missing and duplicated combinations are listed too. A gap is a result: it says a calculation was
attempted and did not come back.

## Figures: the pose worth showing

A trajectory is one pose per disc -- sixty-eight of them for a thirteen angstrom tunnel. Taking
five at even spacing is easy and says nothing, because the pose a figure exists to show is the one
at the top of the energy profile, and even spacing lands on it by luck.

```bash
caver-pymol CaverWEB/8HTB/met3in4ywxjawqzf_results.zip \
    --session CaverWEB/8HTB/pymol_qyj16v/pymol_8HTB_renombrado.pml \
    --tunnel-object tun_cl_3 -o poses.pml
```

Then in PyMOL, with the session already loaded: `@poses.pml`.

The object name is looked up by hash. CaverWeb loads the trajectory as `ligand_<hash>.pdbqt` and
the archive is named after the same hash, so renaming the objects to something readable does not
break the link.

What comes out is a plain script with the states already resolved and the reason beside each:

```
# highest energy on the profile: the obstacle | disc 8, 1.5 A, -0.50 kcal/mol
create snap_2_barrier, MethylEsterT3In, 9, 1
show sticks, snap_2_barrier
color cyan, snap_2_barrier and elem C
pseudoatom pose_label_2, selection=(snap_2_barrier), label="barrier - disc 8, 1.5 A, -0.50 kcal/mol"
```

Editing it needs no Python: to drop a pose, delete its lines; to use a different one, change the
number. State N is profile point N, counting from the mouth of the tunnel -- verified against real
trajectories, where 68 discs give exactly 68 states.

It also marks where the route starts and ends, labels each pose with its disc, distance and energy,
and redraws the tunnel as a mesh, since a solid surface hides the ligand inside it. `--extra N`
adds context poses between the three that matter.

## As a library

```python
from caver_translate import scan, rows

tunnels, jobs = scan("CaverWEB/")
for r in rows(tunnels, jobs):
    print(r["receptor"], r["ligand"], r["tunnel"], r["Ea"], r["flags"])
```

## What this is not

It does not run CAVER or CaverDock — it reads what they produced. And the energies it tabulates
**compare, they do not measure**: CaverDock reports approximate docking energies along a path, not
binding free energies, with a rigid receptor unless flexible side chains were enabled.

## Credit

The analysis it automates is described in *CAVER 3.0* (Chovancová et al., PLoS Comput Biol 2012)
and *CaverDock* (Filipovič et al.; Vávra et al., Bioinformatics 2019). Cite them, not this.

## Licence

GPL-3.0-or-later.
