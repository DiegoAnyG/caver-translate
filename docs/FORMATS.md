# What CaverWeb actually produces

Every statement here was checked against real downloads (8HTB, 4D44, 3SQY — 84 calculations), not
taken from documentation. Where the documentation and the files disagreed, the files won, and the
disagreement is noted.

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
