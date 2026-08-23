"""Choosing poses by the profile instead of by counting.

Sampling a trajectory at even spacing is easy and says nothing: the pose a figure exists to show is
the one at the top of the energy profile, and even spacing lands on it by luck. These fix that the
barrier is always picked, that the state numbers line up with the profile, and that two poses that
would draw the same picture are not both drawn.
"""
from caver_translate.figures import PALETTE, choose_states, script
from caver_translate.parse import Point


def points(triples):
    return [Point(distance=d, disc=i, radius=r, energy_lb=e, energy_ub_min=e, energy_ub_max=e)
            for i, (d, r, e) in enumerate(triples)]


ENTERING = points([(0.0, 1.5, -3.2), (0.8, 1.6, -1.0), (1.6, 1.8, -0.5),
                   (6.0, 2.0, -4.0), (13.0, 2.7, -6.7)])


def tags(states):
    return [t for _s, t, _r in states]


def test_the_barrier_is_always_one_of_them():
    """It is the point of the figure, and an evenly spaced sample would have missed it here."""
    states = choose_states(ENTERING)
    assert "barrier" in tags(states)
    barrier = next(s for s, t, _ in states if t == "barrier")
    assert barrier == 3, "state numbering is one-based: profile point 3 is state 3"


def test_states_are_one_based_like_pymol_counts_them():
    states = choose_states(ENTERING)
    assert min(s for s, _t, _r in states) == 1
    assert max(s for s, _t, _r in states) == len(ENTERING)


def test_leaving_starts_at_the_other_end():
    """An out run stores the profile reversed, so the mouth is the last point, not the first."""
    states = choose_states(list(reversed(ENTERING)), bound="first")
    start = next(s for s, t, _ in states if t == "start")
    end = next(s for s, t, _ in states if t == "end")
    assert start == len(ENTERING) and end == 1


def test_the_reason_travels_with_the_state():
    """The whole point: the script says why each pose is there, so it can be edited without help."""
    for _state, _tag, reason in choose_states(ENTERING):
        assert "disc" in reason and "kcal/mol" in reason


def test_a_minimum_at_the_binding_site_is_not_drawn_twice():
    """The deepest point usually is the bound end; two poses on one spot is a wasted colour."""
    flat_end = points([(0.0, 1.5, -3.0), (5.0, 2.0, -1.0), (10.0, 2.7, -6.1), (13.0, 2.7, -6.1)])
    assert "lowest" not in tags(choose_states(flat_end))


def test_a_minimum_that_is_not_the_binding_site_gets_its_own_pose():
    deep_middle = points([(0.0, 1.5, -3.0), (5.0, 2.0, -9.0), (13.0, 2.7, -4.0)])
    assert "lowest" in tags(choose_states(deep_middle))


def test_context_poses_are_added_only_when_asked():
    assert len(choose_states(ENTERING)) < len(choose_states(ENTERING, extra=2))


def test_an_empty_profile_produces_no_script():
    assert choose_states([]) == []
    assert "nothing to draw" in script("Obj", [])


def test_the_script_creates_one_object_per_state_and_colours_it():
    text = script("MethylEsterT3In", ENTERING, tunnel_obj="tun_cl_3")
    states = choose_states(ENTERING)
    assert text.count("create ") == len(states)
    for i, (state, tag, _r) in enumerate(states):
        assert f"create snap_{i + 1}_{tag}, MethylEsterT3In, {state}, 1" in text
        assert PALETTE[i] in text


def test_the_tunnel_is_drawn_as_a_mesh():
    """A solid surface hides the very thing the figure is of."""
    text = script("Obj", ENTERING, tunnel_obj="tun_cl_3")
    assert "hide surface, tun_cl_3" in text
    assert "show mesh, tun_cl_3" in text


def test_it_cleans_up_after_the_previous_run():
    """Run twice without this and the second figure has the first one still in it."""
    text = script("Obj", ENTERING)
    assert "delete snap_*" in text


def test_start_and_end_markers_are_placed():
    text = script("Obj", ENTERING)
    assert "route_start" in text and "route_end" in text
    assert "show spheres, route_start route_end" in text
