# What CaverWeb and CaverDock actually produce

Every statement here was checked against real files — 84 calculations downloaded from CaverWeb
(8HTB, 4D44, 3SQY), and one of them re-run on this machine three ways — not taken from
documentation. Where the documentation and the files disagreed, the files won, and the
disagreement is noted.

Everything up to "What CaverDock writes when you run it yourself" is about a CaverWeb download;
that last section is about running it here.

## The folder

```
CaverWEB/
  8HTB/                                  one sub-folder per receptor
    qyj16v_summary.txt                   CAVER: the tunnel table
    qyj16v_config.txt                    CAVER: probe radius, starting point, seed
    qyj16v_caver_out.zip                 CAVER: raw output
    benzo3inuh6v4fpkxc_results.zip       one CaverDock job
    ...
    pymol_qyj16v/                        a PyMOL session
      pymol.pml                          CaverWeb's own script
      pymol_8HTB_renombrado.pml          a hand-renamed copy
      data/                              structure, tunnels, trajectories
```

The receptor's name survives only as the folder name. Nothing inside carries it.

## Job identity: the hash is the key

An individually exported job is named `<ligand><tunnel><direction><hash>_results.zip`, for example
`benzo3inuh6v4fpkxc_results.zip` — ligand `benzo`, tunnel 3, direction `in`.

**That same hash names the trajectory in the PyMOL session**: `trajectory/ligand_uh6v4fpkxc.pdbqt`.
Renaming the objects to something readable rewrites the object name and leaves the file name, so
the hash is the one thing that survives and is how the energy profile and the 3D object are joined.

A browser's second download is `..._results (1).zip`. It does **not** end in `_results.zip`, so a
glob for that pattern skips it silently and a whole calculation vanishes from the report. Match by
pattern instead.

Bulk PyMOL downloads reuse identifiers and the name can be wrong. See "direction" below.

## `<hash>_results.zip`

```
results.csv     one line: name, vina_affinity, max_energy, min_energy, min_energy_distance, min_energy_disc
results.json    the full profile -- the only file needed
ligand/
  analysis-lb.pdbqt        poses, lower bound, one MODEL per disc
  analysis-ub.pdbqt        poses, upper bound
  analysis-min.pdbqt       the single lowest-energy pose
  analysis-lb.dat          the profile as plain columns
  analysis-ranges.dat      per-disc min/max/count
  caverdock.config         the parameters used, including the docking box centre
  tunnel.dsd               the discretised tunnel
  caverdock.log.0/.1       the run log
```

`results.json` is a list of entries, one per ligand:

```json
[{"name": "ligand", "affinity": null, "hasUb": true,
  "profile": [{"distance": 0.0, "disc": 0, "radius": 2.7,
               "energyLb": -5.8, "energyUbMin": -5.9, "energyUbMax": -5.9}],
  "maxPoint": {...}, "minPoint": {...}}]
```

`disc` is **not** contiguous — discs can be skipped (0, 1, 2, 3, 5). Index the profile by position,
never by disc number.

A combination that failed produces an archive with no `results.json`, and CaverWeb writes no log
for it. The empty archive and the gap in the table are the only evidence it was ever attempted.

## Which end of the profile is the binding site

**The documentation says disc 0 is the binding site. In the real files it is the surface.** Checked
against the worked example in the project's own interpretation guide and confirmed by the radius,
which is geometry and cannot be mislabelled:

| direction | first disc | last disc |
|---|---|---|
| `in` | radius 1.50 (mouth) | radius 2.70 (cavity) |
| `out` | radius 2.70 (cavity) | radius 1.50 (mouth) |

A tunnel is narrow at its mouth and opens into the cavity, so **the wider end is the binding
site**. The direction reverses the stored profile. `orientation_from_radius` decides from the
radii and flags a disagreement with the file name rather than resolving it silently — that check
found a genuinely mislabelled archive in 3SQY, where two files both named `et1in...` ran in
opposite directions.

## `<hash>_summary.txt` — the CAVER tunnel table

A legend explaining every column, then a fixed-width table under a header beginning `ID  No
No_snaps`. The legend repeats the column names in prose, so the table has to be found by its header
line, not by matching words.

Columns used: `ID`, `Avg_BR` (bottleneck radius), `Avg_L` (length), `Avg_C` (curvature),
`Priority`.

## PyMOL states map one-to-one to profile points

Verified on real trajectories:

| object | profile points | MODEL records |
|---|---|---|
| MethylEsterT3In | 68 | 68 |
| AcidT3In | 68 | 68 |
| AcidT1In | 5 | 5 |

So PyMOL state N (one-based) is profile point N. The states are steps along one route, not
alternatives to choose between.

Atom counts per model identify the compound: acid 17, methyl ester 20, ethyl 23, propyl 26, pentyl
32 — with all atoms. PDBQT merges non-polar hydrogens, so the same compound counts fewer there
(14, 14, 15, 16, 18). Both columns are correct; neither is an error.

## The five reported numbers

```
E_surface   at the mouth
E_max       the barrier, highest anywhere on the profile
E_bound     at the active site

Ea     = E_max   - E_surface     what entering costs
dE_BS  = E_bound - E_surface     how much better the destination is
```

Worked example, 8HTB tunnel 3, benzofuroxanic acid, entering: E_surface −3.2, E_max −0.5,
E_bound −6.7, so Ea = 2.7 and dE_BS = −3.5. That case is a unit test.

## What misleads these numbers

- **A tunnel too short to obstruct anything.** 8HTB tunnel 1: priority 0.946, the best of the six,
  bottleneck radius 2.66 — and 1.00 A long. It is the mouth of the pocket. Priority and length have
  to be read together.
- **A positive E_surface.** The ligand already clashes at the entrance, so dE_BS looks excellent
  because a positive number was subtracted.
- **Mixing lower and upper bound.** LB constrains position only and can pass through
  discontinuities, understating a barrier. UB constrains rotation too. Never compare across them.
- **Assuming in and out are symmetric.** In DHFR tunnel 4 the entering Ea was 2.2–3.9 and the
  leaving Ea reached 10.6–13.4.

## The PyMOL session

`pymol.pml` begins with `cd data`, so PyMOL must be standing in the folder that contains `data/`.
Run from anywhere else and the structure and the tunnel meshes fail to load while the ligands
succeed, leaving a scene that looks empty for no stated reason.

Objects: `structure` (the protein), `tun_cl_1`..`tun_cl_N` (the tunnel clusters), and one object
per trajectory. `set all_states, 0` in the session shows a single state; `set all_states, on` shows
the whole route.

## What CaverDock writes when you run it yourself

`results.json` is CaverWeb's own summary. CaverDock never writes it, so a run made locally has to
be read from the files the engine and its wrappers leave behind. All of it was checked against the
same reference job — 8HTB, benzofuroxanic acid, tunnel 3 — so the two routes can be compared
number for number.

Three layouts turn up, and one reader handles all three because it looks for the files rather than
for the folder shape.

| | what runs it | what it leaves |
|---|---|---|
| a bare run | `mpirun -np 4 caverdock … --out out` | `out-lb.pdbqt`, `out-ub.pdbqt`, `out-ranges.dat` |
| one job | `cd-analysis` | the above plus `analysis-lb.dat`, `analysis-ub.dat`, `tunnel.dsd`, logs |
| a batch | `cd-screening` | `experiments/<name>/{inputs,intermediate,results}`, a shared `cache/` |

### The trajectory

`<name>-lb.pdbqt` is one `MODEL` per disc, in order. Inside each:

```
MODEL 2
REMARK CAVERDOCK RESULT:      -2.7      0.000      0.000
REMARK CAVERDOCK CONSTRAINTS:       0.0      0.0      0.0      0.0
REMARK CAVERDOCK TUNNEL: 1     -2.7      1.6      2.1
```

`TUNNEL` is `<disc> <energy> <radius> <n>`. **The first `TUNNEL` remark of a `MODEL` is that
model's disc.** A run can append a second one repeating disc 0 as a reference, and counting remarks
instead reads that reference as data — 136 remarks for a 68-disc tunnel. Taking the first of each
model reproduces `results.json` exactly for the reference job: all 68 discs, energy and radius
alike.

The fourth column is **not** the distance. It runs 1.3, 2.1, 3.3 … 68.4 over a tunnel 13.0 Å long.

### The energy profile

`.dat` files carry six columns. The names are CaverDock's own, from
`pycaverdock/energy_profile.py`:

```
distance  disc  minE  maxE  radius  lbE
```

For a `cd-analysis` result every column is sound: `minE`, `radius` and the distance match
`results.json` on all 68 discs with nothing left over, and `-ub.dat` supplies `energyUbMin` and
`energyUbMax`.

**In a `cd-screening` `profile.dat` they are not.** `radius` and `lbE` hold one value repeated down
every row, and `minE` comes out as `min(trajectory, lbE)` — 78 of 79 rows. That clips the profile
at the free-docking energy exactly where it rises, which is at the barrier:

| | from `profile.dat` | from the trajectory |
|---|---|---|
| E_surface | −4.6 | −4.6 |
| E_bound | −6.5 | −6.5 |
| **E_max** | **−3.2** | **−2.0** |
| **Ea** | **1.4** | **2.6** |

cd-screening's own `results.csv` is built from those columns, and for that run reported
`E_bound = E_max = E_surface = −3.2` and `dE_BS = 0.0` — a route with no barrier and no binding.
The calculation was fine; the summary of it was not. So energies and radii are read from the
trajectory, and the `.dat` is read for its distance column, which was sound in every file seen.

Only the lower-bound screening has been measured this way. An upper-bound one is read the same way
and has not been checked.

### The upper bound is not in the trajectory

`-ub.pdbqt` is the search, not the result. It revisits discs — 174 remarks over 68 discs, five of
them visited more than once — and neither the first energy recorded for a disc nor the highest
matches the resolved profile: 59 of 68 disagree. Upper-bound numbers come from a `.dat` and from
nowhere else. Without one the job is lower-bound only, and is reported as such.

`cd-energyprofile -d tunnel.dsd -t <name>-ub.pdbqt` writes one, to stdout.

### Distance, when there is no `.dat`

The discretised tunnel `*.dsd` is one row per disc: centre *(x, y, z)*, normal, radius. CaverDock
measures each step **along the disc normal**, not centre to centre — the difference is everywhere
the tunnel bends. Against the reference profile:

| | end to end |
|---|---|
| projected on the normal | 13.0163 Å |
| CaverWeb | 13.0186 Å |
| plain centre-to-centre sum | 13.474 Å |

A run can leave more than one `.dsd`: cd-screening keeps the discretised tunnel and then the
extended one it actually docked into. The one with a disc for every point of the profile is that
one.

### Disc numbers cannot be joined on

Both formats misnumber the last disc. `results.json` ends at 68 for a trajectory that ends at 67;
a cd-screening `profile.dat` skips 78 and ends at 79. The profile and the trajectory are therefore
lined up by position, with one row per model as the check that they correspond at all.

### What names the calculation

Nothing inside the output does. `cd-screening` puts it in the folder name —
`r<receptor>-l<ligand>-t<tunnel>-d<direction>-<lowerbound|upperbound>`, from the basenames of the
three inputs — and that is the only record of what was calculated. A folder named anything else
contributes what it can and leaves the rest unknown rather than inventing it.

### Two things that bite

`cd-screening --no-progress`, the flag its own help recommends for schedulers, crashes on the first
experiment: `AttributeError: 'NoneType' object has no attribute 'update'` in `screening.py:75`,
because the progress bar is used whether or not it was created. Run it without the flag.

One screening block is one direction and one trajectory type. Entering and leaving, lower and upper
bound, is four blocks in the YAML, not one.
