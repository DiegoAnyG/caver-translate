"""The five numbers, and the four ways they mislead.

Every trap tested here was met in a real dataset before it was written down: a tunnel one angstrom
long outscoring every other route, a profile whose mouth energy was positive so the site looked
wonderful by subtraction, and two archives with the same name running in opposite directions.
"""
from caver_translate.metrics import coverage, evaluate, orientation_from_radius
from caver_translate.parse import Job, Point, Tunnel


def points(triples):
    return [Point(distance=d, disc=i, radius=r, energy_lb=e, energy_ub_min=e, energy_ub_max=e)
            for i, (d, r, e) in enumerate(triples)]


def job(triples, direction="in", tunnel=3, has_ub=True):
    return Job(receptor="8HTB", ligand="benzo", tunnel=tunnel, direction=direction,
               source="benzo3in_results.zip", profile=points(triples), has_ub=has_ub)


LONG = Tunnel("8HTB", 3, bottleneck_radius=1.49, length=13.68, curvature=1.23, priority=0.71)
STUB = Tunnel("8HTB", 1, bottleneck_radius=2.66, length=1.00, curvature=1.00, priority=0.95)


def test_the_wide_end_is_the_binding_site():
    """A tunnel is narrow at its mouth and opens into the cavity; the geometry settles it."""
    entering = points([(0.0, 1.5, -3.2), (1.6, 1.8, -0.5), (13.0, 2.7, -6.7)])
    assert orientation_from_radius(entering) == "last"
    assert orientation_from_radius(list(reversed(entering))) == "first"


def test_the_worked_example_from_the_real_data():
    """8HTB, benzofuroxanic acid, tunnel 3, entering: the numbers this has to reproduce."""
    m = evaluate(job([(0.0, 1.5, -3.2), (1.6, 1.8, -0.5), (13.0, 2.7, -6.7)]), LONG)
    assert m.energy_surface == -3.2
    assert m.energy_max == -0.5
    assert m.energy_bound == -6.7
    assert round(m.activation, 2) == 2.7        # Ea = E_max - E_surface
    assert round(m.delta_bs, 2) == -3.5         # dE_BS = E_bound - E_surface
    assert m.flags == ()


def test_leaving_is_read_from_the_same_end_as_entering():
    """The out run stores the profile reversed. Read naively, Ea comes out of the wrong end."""
    leaving = job([(0.0, 2.7, -6.8), (11.0, 1.8, -0.4), (13.0, 1.5, -1.9)], direction="out")
    m = evaluate(leaving, LONG)
    assert m.energy_surface == -1.9 and m.energy_bound == -6.8
    assert round(m.activation, 2) == 1.5
    assert "direction_mismatch" not in m.flags


def test_a_name_that_disagrees_with_the_geometry_is_flagged():
    """Two archives in the real data both claimed 'in'; their radii ran opposite ways."""
    mislabelled = job([(0.0, 3.0, -5.9), (1.0, 2.8, -5.1)], direction="in")
    m = evaluate(mislabelled, LONG)
    assert "direction_mismatch" in m.flags
    assert m.orientation == "first"             # the radius was believed, not the file name


def test_a_tunnel_with_nothing_to_cross_is_flagged():
    """CAVER gives it the best priority because there is no distance in which to be obstructed."""
    m = evaluate(job([(0.0, 2.7, -5.8), (1.0, 2.7, -6.4)], tunnel=1), STUB)
    assert "short_tunnel" in m.flags


def test_a_positive_mouth_energy_is_flagged():
    """dE_BS then looks excellent because a positive number was subtracted, not because it binds."""
    m = evaluate(job([(0.0, 1.2, 3.1), (5.0, 1.5, 3.4), (10.0, 2.7, -5.4)]), LONG)
    assert "positive_surface" in m.flags
    assert m.delta_bs < -8                      # the misleading figure, reported with its warning


def test_a_lower_bound_only_run_says_so():
    m = evaluate(job([(0.0, 1.5, -3.0), (5.0, 2.7, -6.0)], has_ub=False), LONG)
    assert "lower_bound_only" in m.flags


def test_a_constant_radius_falls_back_to_the_name_and_admits_it():
    m = evaluate(job([(0.0, 2.0, -3.0), (5.0, 2.0, -6.0)], direction="in"), LONG)
    assert "orientation_from_name" in m.flags
    assert m.energy_bound == -6.0


def test_an_empty_profile_is_a_failure_not_a_zero():
    m = evaluate(Job("8HTB", "pent", 2, "out", "pent2out_results.zip"))
    assert m.flags == ("failed",)
    assert m.activation is None


def test_coverage_names_what_never_came_back():
    """CaverWeb writes no log for a failed combination; the gap is the whole evidence."""
    jobs = [job([(0.0, 1.0, -1.0)], direction="in", tunnel=1),
            job([(0.0, 1.0, -1.0)], direction="out", tunnel=1),
            job([(0.0, 1.0, -1.0)], direction="in", tunnel=2)]
    cov = coverage(jobs)
    assert cov["expected"] == 4 and cov["present"] == 3
    assert ("8HTB", "benzo", 2, "out") in cov["missing"]


def test_coverage_notices_two_archives_claiming_one_combination():
    a = job([(0.0, 1.0, -1.0)], direction="in", tunnel=1)
    b = job([(0.0, 1.0, -1.0)], direction="in", tunnel=1)
    b.source = "et1inother_results.zip"
    cov = coverage([a, b])
    assert cov["duplicated"], "one of the two identifiers was reused and both cannot be right"
