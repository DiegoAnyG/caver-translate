"""An upper bound that was asked for and refused is not the same as one never asked for.

CaverDock writes `<name>-failed.pdbqt` when a trajectory cannot be completed: one model, at the
disc it could not get past, with the log repeating `TRAJECTORY: convergence ... FAILED (for
backtracking N)`. The ligand does not cross there with its rotation constrained.

Read as an absence, that becomes `lower_bound_only` -- the same thing a run says when nobody asked
for an upper bound. One is a question not put, the other is an answer of "it does not pass", and
the second is the more interesting of the two.
"""
from caver_translate.local import failed_at, parse_local_job
from caver_translate.metrics import evaluate

DISCS = [(0, -3.2, 1.5), (1, -2.7, 1.6), (2, -0.5, 1.7), (3, -6.7, 2.7)]
DSD = "\n".join(f"{i}.0 0.0 0.0 1.0 0.0 0.0 {1.5 + i * 0.1}" for i in range(4)) + "\n"


def model(disc, energy, radius):
    return "\n".join([
        f"MODEL {disc + 1}",
        f"REMARK CAVERDOCK TUNNEL: {disc}     {energy}      {radius}      1.3",
        "ENDMDL"])


def run(folder, failed_disc=None):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "analysis-lb.pdbqt").write_text(
        "\n".join(model(d, e, r) for d, e, r in DISCS) + "\n")
    (folder / "tunnel.dsd").write_text(DSD)
    if failed_disc is not None:
        (folder / "analysis-failed.pdbqt").write_text(model(failed_disc, -6.9, 2.8) + "\n")
    return folder


def test_the_disc_it_stopped_at_is_read(tmp_path):
    assert failed_at(run(tmp_path / "a", failed_disc=80)) == 80


def test_a_finished_run_has_no_such_marker(tmp_path):
    assert failed_at(run(tmp_path / "b")) is None


def test_a_refused_upper_bound_is_flagged_as_refused(tmp_path):
    job = parse_local_job(run(tmp_path / "c", failed_disc=80))
    assert job.ub_failed_at == 80
    assert "disc 80" in job.note
    flags = evaluate(job).flags
    assert "upper_bound_failed" in flags
    assert "lower_bound_only" not in flags, "a refusal was reported as a gap"


def test_never_asking_still_reads_as_never_asked(tmp_path):
    job = parse_local_job(run(tmp_path / "d"))
    assert job.ub_failed_at is None
    assert "lower_bound_only" in evaluate(job).flags


def test_the_flag_has_something_to_say_for_itself():
    from caver_translate.report import FLAG_TEXT

    assert "upper_bound_failed" in FLAG_TEXT
    assert "did not converge" in FLAG_TEXT["upper_bound_failed"]
