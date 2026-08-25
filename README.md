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

### Which route a compound takes

The ranking above answers "what is the easiest thing in this dataset". The question actually being
asked is usually narrower: *this* compound, which tunnel does it use, and *this* tunnel, which
compounds get through it. The page carries both, as the same routes grouped two ways.

A route is one receptor, one compound and one tunnel, with entering and leaving side by side --
they are two calculations of the same route, and split across a ranking the same tunnel turns up
twice and gets compared with itself.

How to read a row:

- **`Ea`** is what entering costs. Near zero is a way in; a high one is a real obstacle. It is the
  metric that compares tunnels, because it is the one that says which offers least resistance.
- **`dE_BS`** is how much better the site is than the surface. Very negative means the ligand is
  far more likely inside than out, so binding is stronger in the site.
- **`E_surface`** is the energy at the mouth, and it is what decides whether `dE_BS` can be read at
  all. A positive mouth means the ligand already clashes on arrival, and subtracting a positive
  number makes `dE_BS` look excellent when nothing favourable was measured. Rows like that are
  marked `positive_surface` and their `dE_BS` is printed as *mouth clashes*.
- **Length** belongs beside `Ea`. A tunnel with nothing to cross costs nothing to cross, which is
  why a one angstrom "tunnel" tops a ranking without meaning anything. Those are marked
  `short_tunnel`.

No verdict column, and no score combining the four numbers into one. Where the cut between "enters
easily" and "has an obstacle" falls is a decision about the chemistry, and it is not the tool's to
make: the ranking, the numbers and the marks are on the page, and the reading is yours.


## Figures: the pose worth showing

A trajectory is one pose per disc -- sixty-eight of them for a thirteen angstrom tunnel. Taking
five at even spacing is easy and says nothing, because the pose a figure exists to show is the one
at the top of the energy profile, and even spacing lands on it by luck.

```bash
caver-pymol CaverWEB/8HTB/met3in4ywxjawqzf_results.zip \
    --session CaverWEB/8HTB/pymol_qyj16v/pymol_8HTB_renombrado.pml -o poses.pml
```

### Choosing which ones

Run it with **no arguments** and it asks. Everything on the menus is read off the disk, so there is
nothing to configure and nothing that can go stale:

```
$ cd CaverWEB
$ caver-pymol

Results folders:
    1   3SQY
    2   4D44
    3   8HTB
  > 3

PyMOL session to take the object names from:
    1   pymol_qyj16v/pymol.pml
    2   pymol_qyj16v/pymol_8HTB_renombrado.pml
  > 2

Tunnel:
    1   tunnel 1   (10 trajectories)
    2   tunnel 2   (10 trajectories)
    3   tunnel 3   (10 trajectories)
  > (Enter for all) 3

Compound:
    1   benzo   (6 trajectories)
    2   et      (6 trajectories)
  ...
  > (Enter for all) 1

Direction:
    i   in    -- towards the active site
    o   out   -- away from it
  > (Enter for all) i

1 of 30 trajectories. The same run without the questions:
  caver-pymol benzo3inuh6v4fpkxc_results.zip --session .../pymol_8HTB_renombrado.pml -o .../poses
benzo3inuh6v4fpkxc_results.zip  ->  @C:/.../8HTB/poses/benzo3inuh6v4fpkxc.pml
```

Enter on any question takes everything it lists, so tunnel 3 and Enter for the compound is the same
tunnel across all five, in one run. The last line is the command that would have done it without
the questions, printed so that the second time can skip them.

### Choosing them without the questions

The archive names already say the compound, the tunnel and the direction
(`benzo3inuh6v4fpkxc_results.zip` is the acid, tunnel 3, entering), so the shell already knows how
to select them and there is no set of flags to learn. Name as many archives as you like, or a
folder, and give `-o` a folder to write into:

```bash
cd CaverWEB/8HTB
S=pymol_qyj16v/pymol_8HTB_renombrado.pml

caver-pymol *3in*_results*.zip   --session $S -o poses/   # one tunnel, every compound
caver-pymol benzo*_results*.zip  --session $S -o poses/   # one compound, every tunnel
caver-pymol *out*_results*.zip   --session $S -o poses/   # everything leaving
caver-pymol .                    --session $S -o poses/   # the lot
```

One `.pml` per archive, named after it. The tunnel is read from the archive name, so
`--tunnel-object` is only needed to override it (`--tunnel-object ''` draws no tunnel at all). The
object name is looked up by hash, which is what makes a whole folder possible: none of the thirty
names has to be typed.

### Using it

The script draws poses **into a session that is already open**. On its own it has nothing to draw
on, so load the CaverWeb session first, from the folder that holds `data/`:

```
cd C:/.../CaverWEB/8HTB/pymol_qyj16v
@pymol_8HTB_renombrado.pml
@C:/.../CaverWEB/8HTB/pymol_qyj16v/poses/met3in4ywxjawqzf.pml
```

That first `cd` matters: the session script starts with `cd data`, so PyMOL has to be standing in
the folder above it. Run it from anywhere else and the structure and the tunnels never load.

The full path on the third line matters for the same reason, backwards: the session script ends
inside `data/` and never comes back, so a relative `@poses/...` is looked for in `data/poses/` and
is not found. `caver-pymol` prints the absolute path for each script it writes; paste that.

Load the session **once**. PyMOL appends states, it does not replace them, so a second
`@pymol_....pml` leaves every object holding two copies of itself: memory runs out, undo switches
itself off, and the object panel starts throwing errors that name no cause. The generated script
counts the states it finds against the profile and says so if this has happened.

If the session is not there, the script says so once and names what is missing, instead of failing
on every line.

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

It also clears the way so the route can be seen: only the tunnel this trajectory goes through,
drawn as a mesh rather than a solid surface; the receptor faded to 80 % transparent with its waters
removed; the other tunnels switched off but still loaded (`enable tun_cl_*` brings them back).

Each pose carries the two numbers a figure is read for -- how far it is from the active site and
what the energy is there -- and the barrier says what it is:

```
entrance       13.0 A from site   -5.4 kcal/mol
HARDEST STEP   11.0 A from site   -3.2 kcal/mol
deepest point   4.5 A from site   -6.9 kcal/mol
active site     0.0 A from site   -6.4 kcal/mol
```

`--extra N` adds context poses between them; `--receptor-object ""` leaves the protein alone.

## Runs made on your own machine

CaverDock does not write `results.json` — that is CaverWeb's own summary. Point it at what
CaverDock did write and it reads that instead, whether the job was a bare `caverdock` run, one
`cd-analysis`, or a `cd-screening` batch:

```bash
caver-translate screening_out/ -o report/
```

Two things it does rather than take the files as written, both because a real run showed why:

- **Energies and radii come from the trajectory**, not from the profile `.dat`. They agree on a
  `cd-analysis` result — all 68 discs of the reference job, exactly — but a `cd-screening`
  `profile.dat` clips its lower bound at the free-docking energy, which flattens the barrier and
  nothing else. On the run measured here that moved `E_max` from −2.0 to −3.2 and `Ea` from 2.6 to
  1.4, and cd-screening's own `results.csv`, built from the same columns, reported a route with no
  barrier and no binding. The calculation was fine; the summary of it was not.
- **Distances are measured along the disc normal** when there is no `.dat` to take them from, which
  is how CaverDock measures them. The plain centre-to-centre sum runs 3.5 % long.

Upper bounds come from a `-ub.dat` and nowhere else. The upper-bound trajectory is the search
rather than the result, and 59 of 68 discs disagree with the resolved profile; without the file the
job is reported `lower_bound_only`. `cd-energyprofile -d tunnel.dsd -t out-ub.pdbqt` writes one.

## As a library

```python
from caver_translate import scan, rows

tunnels, jobs = scan("CaverWEB/")
for r in rows(tunnels, jobs):
    print(r["receptor"], r["ligand"], r["tunnel"], r["Ea"], r["flags"])
```

## The file formats

Everything verified about what CaverWeb and CaverDock produce -- the folder, the archives, the
files a local run leaves, which end of the profile is the binding site, and the four things that
mislead the numbers -- is in
[docs/FORMATS.md](docs/FORMATS.md). It was checked against real downloads, and where the upstream
documentation disagreed with the files, the files won.

## What this is not

It does not run CAVER or CaverDock — it reads what they produced. And the energies it tabulates
**compare, they do not measure**: CaverDock reports approximate docking energies along a path, not
binding free energies, with a rigid receptor unless flexible side chains were enabled.

## Credit

The analysis it automates is described in *CAVER 3.0* (Chovancová et al., PLoS Comput Biol 2012)
and *CaverDock* (Filipovič et al.; Vávra et al., Bioinformatics 2019). Cite them, not this.

## Licence

GPL-3.0-or-later.
