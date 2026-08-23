"""Reading the two formats CaverWeb hands back.

Fixtures are built here rather than checked in: the real downloads are megabytes of poses and
logs, and every fact these tests rely on is a few lines of the two files that matter.
"""
import json
import zipfile
from pathlib import Path

from caver_translate.parse import parse_job, parse_tunnels, scan

# The CAVER summary is a legend followed by a fixed-width table. The legend repeats the column
# names in prose, which is why the table is found by its header line and not by matching words.
SUMMARY = """            ID: Identification of a given tunnel cluster = ranking of a given cluster
        Avg_BR: Average bottleneck radius [A].
         Avg_L: Average tunnel length [A].
      Priority: Tunnel priority calculated by averaging tunnel throughputs.

  ID      No   No_snaps   Avg_BR       SD   Max_BR    Avg_L      SD   Avg_C      SD    Priority
   1       1          1    2.661    0.000     2.66    1.000   0.000   1.000   0.000     0.94646
   2       1          1    1.522    0.000     1.52    8.814   0.000   1.335   0.000     0.74064


-----------------------------------------------------------------------
 Thank you for using CAVER, please cite:
"""


def profile(points):
    return [{"distance": d, "disc": i, "radius": r, "energyLb": e,
             "energyUbMin": e, "energyUbMax": e} for i, (d, r, e) in enumerate(points)]


def make_zip(folder, name, points, has_ub=True, entry_name="ligand"):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    body = [{"name": entry_name, "affinity": None, "hasUb": has_ub, "profile": profile(points)}]
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("results.json", json.dumps(body))
        zf.writestr("ligand/analysis-lb.dat", "")
    return path


def test_the_tunnel_table_is_read_and_the_legend_is_not(tmp_path):
    """The legend names every column in prose; parsing it as data invents tunnels."""
    receptor = tmp_path / "8HTB"
    receptor.mkdir()
    (receptor / "abc_summary.txt").write_text(SUMMARY)
    tunnels = parse_tunnels(receptor / "abc_summary.txt")
    assert len(tunnels) == 2
    first = tunnels[0]
    assert (first.receptor, first.tunnel) == ("8HTB", 1)
    assert first.bottleneck_radius == 2.661
    assert first.length == 1.0
    assert first.priority == 0.94646


def test_the_job_name_gives_ligand_tunnel_and_direction(tmp_path):
    path = make_zip(tmp_path / "8HTB", "benzo3inuh6v4fpkxc_results.zip",
                    [(0.0, 1.5, -3.2), (1.0, 1.8, -0.5), (2.0, 2.7, -6.7)])
    job = parse_job(path)[0]
    assert (job.ligand, job.tunnel, job.direction) == ("benzo", 3, "in")
    assert len(job.profile) == 3
    assert job.has_ub


def test_a_browser_duplicated_download_is_still_understood(tmp_path):
    """Downloading the same job twice gives "name (1).zip"; the job is unchanged."""
    path = make_zip(tmp_path / "3SQY", "benzo1inlvikai9svq_results (1).zip", [(0.0, 1.0, -1.0)])
    job = parse_job(path)[0]
    assert (job.ligand, job.tunnel, job.direction) == ("benzo", 1, "in")


def test_an_archive_without_a_profile_is_reported_not_dropped(tmp_path):
    """A failed CaverWeb combination leaves no log, so the empty archive is the only evidence."""
    folder = tmp_path / "8HTB"
    folder.mkdir(parents=True)
    path = folder / "met2outzzz_results.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("ligand/log.txt", "crashed")
    job = parse_job(path)[0]
    assert job.profile == []
    assert "no results.json" in job.note
    assert job.tunnel == 2


def test_a_corrupt_archive_does_not_stop_the_run(tmp_path):
    folder = tmp_path / "8HTB"
    folder.mkdir(parents=True)
    path = folder / "et1inaaa_results.zip"
    path.write_bytes(b"not a zip at all")
    job = parse_job(path)[0]
    assert "unreadable" in job.note


def test_scan_finds_receptors_as_sub_folders(tmp_path):
    for receptor in ("8HTB", "4D44"):
        (tmp_path / receptor).mkdir()
        (tmp_path / receptor / "x_summary.txt").write_text(SUMMARY)
        make_zip(tmp_path / receptor, "benzo1inaaa_results.zip", [(0.0, 1.0, -1.0), (1.0, 2.0, -2.0)])
    tunnels, jobs = scan(tmp_path)
    assert {t.receptor for t in tunnels} == {"8HTB", "4D44"}
    assert {j.receptor for j in jobs} == {"8HTB", "4D44"}


def test_scan_also_accepts_a_single_receptor_folder(tmp_path):
    """Someone who downloaded one target points at that folder, not at a parent of one."""
    (tmp_path / "x_summary.txt").write_text(SUMMARY)
    make_zip(tmp_path, "benzo1inaaa_results.zip", [(0.0, 1.0, -1.0), (1.0, 2.0, -2.0)])
    tunnels, jobs = scan(tmp_path)
    assert len(tunnels) == 2 and len(jobs) == 1


def test_a_second_download_is_still_found(tmp_path):
    """It is named "..._results (1).zip". Globbing for "*_results.zip" skipped it silently, and a
    whole calculation disappeared from the report -- the one failure this tool must not have."""
    receptor = tmp_path / "3SQY"
    make_zip(receptor, "benzo1inlvikai9svq_results (1).zip", [(0.0, 1.0, -1.0), (1.0, 2.0, -2.0)])
    _tunnels, jobs = scan(tmp_path)
    assert len(jobs) == 1
    assert (jobs[0].ligand, jobs[0].tunnel, jobs[0].direction) == ("benzo", 1, "in")


def test_the_other_archives_in_the_folder_are_left_alone(tmp_path):
    """A CaverWeb folder also holds the CAVER output and a PyMOL session; neither is a job."""
    receptor = tmp_path / "4D44"
    make_zip(receptor, "benzo1inaaa_results.zip", [(0.0, 1.0, -1.0), (1.0, 2.0, -2.0)])
    (receptor / "z9uqie_caver_out.zip").write_bytes(b"PK")
    (receptor / "pymol_z9uqie.zip").write_bytes(b"PK")
    _tunnels, jobs = scan(tmp_path)
    assert len(jobs) == 1
