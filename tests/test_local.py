"""Reading a CaverDock run made here rather than a CaverWeb download.

The numbers in these fixtures are the first discs of the reference job -- 8HTB, benzofuroxanic
acid, tunnel 3 -- so a change that breaks the format shows up as a wrong energy and not just a
missing file.
"""
import math

from caver_translate.local import disc_distances, local_jobs, parse_local_job, trajectory
from caver_translate.metrics import evaluate
from caver_translate.parse import scan

# distance disc energy_min energy_max radius energy_reference
LB_DAT = """0.0 0 -3.2 -3.2 1.5 -3.2
0.20950621109268042 1 -2.7 -2.7 1.6 -2.7
0.4188202511301231 2 -2.3 -2.3 1.7 -2.3
0.6028219537413673 3 -0.5 -0.5 1.7 -0.5
0.8423452955582853 4 -6.7 -6.7 2.7 -6.7
"""

# The same discs of the upper bound. Disc 2 is where min and max differ, which is the whole reason
# the upper bound is read from here and not from the trajectory.
UB_DAT = """0.0 0 -2.9 -2.9 1.5 -3.2
0.20950621109268042 1 -2.7 -2.7 1.6 -2.7
0.4188202511301231 2 -2.5 -0.4 1.7 -2.5
0.6028219537413673 3 0.1 0.1 1.7 -2.3
0.8423452955582853 4 -5.5 -5.5 2.7 -6.7
"""


def model(disc, energy, radius, extra_disc_zero=False):
    """One MODEL of a trajectory, optionally with the trailing disc-0 remark a bare run adds."""
    lines = [f"MODEL {disc + 1}",
             f"REMARK CAVERDOCK RESULT:      {energy}      0.000      0.000",
             "REMARK CAVERDOCK CONSTRAINTS:       0.0      0.0      0.0      0.0",
             f"REMARK CAVERDOCK TUNNEL: {disc}     {energy}      {radius}      1.3"]
    if extra_disc_zero:
        lines += ["REMARK CAVERDOCK RESULT:      -3.2      0.000      0.000",
                  "REMARK CAVERDOCK CONSTRAINTS:       0.0      0.0      0.0      0.0",
                  "REMARK CAVERDOCK TUNNEL: 0     -3.2      1.5      1.3"]
    return "\n".join(lines + ["ENDMDL"])


DISCS = [(0, -3.2, 1.5), (1, -2.7, 1.6), (2, -2.3, 1.7), (3, -0.5, 1.7), (4, -6.7, 2.7)]

# x y z, then the disc normal, then the radius. A straight run along x one angstrom apart.
DSD = "\n".join(f"{i}.0 0.0 0.0 1.0 0.0 0.0 {1.5 + i * 0.1}" for i in range(5)) + "\n"


def write_run(folder, extra_disc_zero=False, dat=False, ub_dat=False, name="out"):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{name}-lb.pdbqt").write_text(
        "\n".join(model(d, e, r, extra_disc_zero) for d, e, r in DISCS) + "\n")
    (folder / "tunnel.dsd").write_text(DSD)
    if dat:
        (folder / f"{name}-lb.dat").write_text(LB_DAT)
    if ub_dat:
        (folder / f"{name}-ub.dat").write_text(UB_DAT)
    return folder


def test_the_first_remark_of_a_model_is_that_models_disc(tmp_path):
    """A bare caverdock run repeats disc 0 inside every MODEL. Counting remarks reads it as data."""
    folder = write_run(tmp_path / "job", extra_disc_zero=True)
    steps = trajectory(folder / "out-lb.pdbqt")
    assert [d for d, _e, _r in steps] == [0, 1, 2, 3, 4]
    assert [e for _d, e, _r in steps] == [-3.2, -2.7, -2.3, -0.5, -6.7]


def test_distance_follows_the_disc_normal(tmp_path):
    """Along a straight tunnel the projection and the centre-to-centre distance agree."""
    folder = write_run(tmp_path / "job")
    assert disc_distances(folder / "tunnel.dsd") == [0.0, 1.0, 2.0, 3.0, 4.0]


def test_a_bare_run_is_read_from_the_trajectory(tmp_path):
    """No .dat anywhere: the profile still has to come out, with the distances measured."""
    job = parse_local_job(write_run(tmp_path / "job", extra_disc_zero=True))
    assert [p.energy_lb for p in job.profile] == [-3.2, -2.7, -2.3, -0.5, -6.7]
    assert [p.radius for p in job.profile] == [1.5, 1.6, 1.7, 1.7, 2.7]
    assert job.profile[-1].distance == 4.0
    assert job.has_ub is False


def test_the_upper_bound_comes_only_from_its_dat_file(tmp_path):
    """The trajectory holds the search, not the resolved bound. Without the file, say so."""
    without = parse_local_job(write_run(tmp_path / "a", dat=True))
    assert without.has_ub is False
    assert "cd-energyprofile" in without.note

    with_ub = parse_local_job(write_run(tmp_path / "b", dat=True, ub_dat=True))
    assert with_ub.has_ub is True
    assert [p.energy_ub_min for p in with_ub.profile] == [-2.9, -2.7, -2.5, 0.1, -5.5]
    assert [p.energy_ub_max for p in with_ub.profile] == [-2.9, -2.7, -0.4, 0.1, -5.5]
    # the lower bound is still the lower bound: -ub.dat must not overwrite it
    assert [p.energy_lb for p in with_ub.profile] == [-3.2, -2.7, -2.3, -0.5, -6.7]


def test_a_missing_dsd_is_reported_not_guessed(tmp_path):
    folder = write_run(tmp_path / "job")
    (folder / "tunnel.dsd").unlink()
    job = parse_local_job(folder)
    assert len(job.profile) == 5
    assert "distances are unknown" in job.note


def test_a_screening_profile_is_read_and_its_radii_repaired(tmp_path):
    """cd-screening writes results/profile.dat with the radius repeated down the column.

    Left as written it reports a tunnel of constant width, which decides the wrong end is the
    binding site. The discs it was docked into are on disk; the radii come from there.
    """
    exp = tmp_path / "experiments" / "r8HTB-lbenzo-ttun_cl_003-dout-lowerbound"
    (exp / "results").mkdir(parents=True)
    # distance disc minE maxE radius lbE -- only minE varies
    (exp / "results" / "profile.dat").write_text(
        "".join(f"{i}.0 {i} {e} -3.2 1.5 -3.2\n" for i, e in enumerate([-6.5, -6.0, -5.0, -3.2, -4.6])))
    # docked outwards, so the discs run from the cavity to the mouth and narrow as they go
    (exp / "intermediate").mkdir()
    (exp / "intermediate" / "tunnel.dsd").write_text(
        "\n".join(f"{i}.0 0.0 0.0 1.0 0.0 0.0 {1.9 - i * 0.1}" for i in range(5)) + "\n")

    job = parse_local_job(exp)
    assert job.tunnel == 3 and job.direction == "out"
    assert [p.energy_lb for p in job.profile] == [-6.5, -6.0, -5.0, -3.2, -4.6]
    assert [p.radius for p in job.profile] == [1.9, 1.8, 1.7, 1.6, 1.5]

    m = evaluate(job)
    assert m.orientation == "first"          # docked outwards: disc 0 is the site
    assert m.energy_bound == -6.5 and m.energy_surface == -4.6
    assert "orientation_from_name" not in m.flags


def test_the_barrier_comes_from_the_trajectory_not_the_profile_file(tmp_path):
    """cd-screening clips minE at the free-docking energy, which flattens the barrier.

    The trajectory beside it is not clipped, and one row per model is what says they line up.
    Reading the profile file instead moved E_max from -0.5 to -3.2 on the run this is taken from.
    """
    exp = tmp_path / "experiments" / "r8HTB-lbenzo-t3-din-lowerbound"
    write_run(exp / "intermediate" / "caverdock", name="caverdock")
    (exp / "results").mkdir(parents=True)
    clipped = [min(e, -3.2) for _d, e, _r in DISCS]
    # the last disc is misnumbered, as both CaverWeb and cd-screening do at the tail
    numbers = [0, 1, 2, 3, 5]
    (exp / "results" / "profile.dat").write_text(
        "".join(f"{i * 2}.0 {n} {e} -3.2 1.5 -3.2\n" for i, (n, e) in enumerate(zip(numbers, clipped))))

    job = parse_local_job(exp)
    assert job.source == "caverdock-lb.pdbqt"
    assert [p.energy_lb for p in job.profile] == [-3.2, -2.7, -2.3, -0.5, -6.7]   # not clipped
    assert [p.distance for p in job.profile] == [0.0, 2.0, 4.0, 6.0, 8.0]         # from the file
    assert evaluate(job).energy_max == -0.5


def test_a_screening_experiment_is_found_by_its_profile(tmp_path):
    exp = tmp_path / "out" / "lb_out" / "experiments" / "r8HTB-lbenzo-t3-dout-lowerbound"
    (exp / "results").mkdir(parents=True)
    (exp / "results" / "profile.dat").write_text("0.0 0 -6.5 -3.2 1.5 -3.2\n1.0 1 -4.6 -3.2 1.5 -3.2\n")
    jobs = local_jobs(tmp_path)
    assert [j.ligand for j in jobs] == ["benzo"]


def test_a_screening_folder_name_says_what_was_calculated(tmp_path):
    """cd-screening records receptor, ligand, tunnel and direction only in the folder name."""
    name = "r8HTB-lbenzo-ttun_cl_003-din-lowerbound"
    job = parse_local_job(write_run(tmp_path / "experiments" / name, dat=True))
    assert (job.receptor, job.ligand, job.tunnel, job.direction) == ("8HTB", "benzo", 3, "in")


def test_an_unnamed_folder_claims_nothing(tmp_path):
    """Better an unknown tunnel than an invented one: coverage counts on it."""
    job = parse_local_job(write_run(tmp_path / "8HTB" / "whatever", dat=True))
    assert (job.receptor, job.ligand) == ("8HTB", "whatever")
    assert job.tunnel is None and job.direction is None


def test_a_tree_of_runs_is_found_whatever_the_nesting(tmp_path):
    write_run(tmp_path / "one", dat=True)
    write_run(tmp_path / "deep" / "deeper" / "two", dat=True)
    assert len(local_jobs(tmp_path)) == 2


def test_scan_falls_through_to_a_local_run(tmp_path):
    """A folder with no archives in it is not empty -- it may be a run made here."""
    write_run(tmp_path / "job", dat=True, ub_dat=True)
    _tunnels, jobs = scan(tmp_path)
    assert len(jobs) == 1
    assert jobs[0].has_ub is True


def test_the_reported_numbers_come_out_of_a_local_run(tmp_path):
    """End to end: the five numbers, from files CaverDock wrote, with no results.json involved."""
    job = parse_local_job(write_run(tmp_path / "job", dat=True, ub_dat=True))
    m = evaluate(job)
    assert m.orientation == "last"          # the tunnel widens towards the site
    assert m.energy_surface == -3.2 and m.energy_bound == -6.7
    assert math.isclose(m.activation, 2.7, abs_tol=1e-9)
    assert math.isclose(m.delta_bs, -3.5, abs_tol=1e-9)
    assert "lower_bound_only" not in m.flags
