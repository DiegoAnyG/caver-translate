"""Choosing poses by the profile instead of by counting, and drawing them so they can be seen.

Sampling a trajectory at even spacing is easy and says nothing: the pose a figure exists to show is
the one at the top of the energy profile, and even spacing lands on it by luck. The rest is getting
everything else out of its way -- the protein, the other tunnels, the markers.
"""
from caver_translate.figures import PALETTE, choose_states, script
from caver_translate.parse import Point


def points(triples):
    return [Point(distance=d, disc=i, radius=r, energy_lb=e, energy_ub_min=e, energy_ub_max=e)
            for i, (d, r, e) in enumerate(triples)]


ENTERING = points([(0.0, 1.5, -3.2), (0.8, 1.6, -1.0), (1.6, 1.8, -0.5),
                   (6.0, 2.0, -4.0), (13.0, 2.7, -6.7)])


def tags(states):
    return [tag for _state, tag, _label, _reason in states]


def labels(states):
    return [label for _state, _tag, label, _reason in states]


def test_the_barrier_is_always_one_of_them():
    """It is the point of the figure, and an evenly spaced sample would have missed it here."""
    states = choose_states(ENTERING)
    assert "barrier" in tags(states)
    barrier = next(s for s, t, _l, _r in states if t == "barrier")
    assert barrier == 3, "state numbering is one-based: profile point 3 is state 3"


def test_states_are_one_based_like_pymol_counts_them():
    states = choose_states(ENTERING)
    assert min(s for s, _t, _l, _r in states) == 1
    assert max(s for s, _t, _l, _r in states) == len(ENTERING)


def test_leaving_starts_at_the_other_end():
    """An out run stores the profile reversed, so the mouth is the last point, not the first."""
    states = choose_states(list(reversed(ENTERING)), bound="first")
    start = next(s for s, t, _l, _r in states if t == "start")
    end = next(s for s, t, _l, _r in states if t == "end")
    assert start == len(ENTERING) and end == 1


def test_distance_is_measured_from_the_active_site():
    """That is the axis the figure is read along. Measured from the file's first point instead, a
    run stored backwards would count from the wrong end and every label would be inverted."""
    entering = choose_states(ENTERING)
    at_site = next(l for _s, t, l, _r in entering if t == "end")
    at_mouth = next(l for _s, t, l, _r in entering if t == "start")
    assert at_site.startswith("active site   0.0 A from site")
    assert "13.0 A from site" in at_mouth

    leaving = choose_states(list(reversed(ENTERING)), bound="first")
    assert next(l for _s, t, l, _r in leaving if t == "end").startswith("active site   0.0 A")
    assert "13.0 A from site" in next(l for _s, t, l, _r in leaving if t == "start")


def test_the_label_carries_the_distance_and_the_energy_and_nothing_else():
    for label in labels(choose_states(ENTERING)):
        assert "A from site" in label and "kcal/mol" in label
        assert "disc" not in label, "the disc number means nothing to a reader of the figure"


def test_the_hardest_step_says_so_on_the_figure():
    barrier = next(l for _s, t, l, _r in choose_states(ENTERING) if t == "barrier")
    assert barrier.startswith("HARDEST STEP")


def test_the_reason_travels_with_the_state():
    """The whole point: the script says why each pose is there, so it can be edited without help."""
    for _state, _tag, _label, reason in choose_states(ENTERING):
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
    for i, (state, tag, _label, _reason) in enumerate(states):
        assert f"create snap_{i + 1}_{tag}, MethylEsterT3In, {state}, 1" in text
        assert PALETTE[i] in text


def test_only_the_tunnel_in_question_is_shown():
    """The session loads every cluster. Six tunnels over one route is not a figure of anything."""
    text = script("Obj", ENTERING, tunnel_obj="tun_cl_3")
    assert "hide everything, tun_cl_*" in text
    assert "enable tun_cl_3" in text
    assert "show mesh, tun_cl_3" in text
    assert "set mesh_width, 0.3" in text


def test_the_receptor_is_faded_and_its_waters_removed():
    """A solid protein hides the tunnel and every pose in it, and each water is another dot."""
    text = script("Obj", ENTERING, receptor_obj="structure")
    assert "set cartoon_transparency, 0.8, structure" in text
    assert "remove structure and solvent" in text


def test_the_receptor_can_be_left_alone():
    assert "cartoon_transparency" not in script("Obj", ENTERING, receptor_obj="")


def test_the_ends_are_marked_with_text_not_with_spheres():
    """A sphere at each end is one more object in front of the route, and the labels already say
    which end is which."""
    text = script("Obj", ENTERING)
    assert "show spheres" not in text
    assert "entrance" in text and "active site" in text


def test_it_cleans_up_after_the_previous_run():
    """Run twice without this and the second figure has the first one still in it."""
    text = script("Obj", ENTERING)
    assert "delete snap_*" in text


def test_the_script_is_seven_bit_ascii():
    """PyMOL read a middle dot in a label as latin-1 and printed mojibake. Same lesson as any
    other file handed to a tool whose encoding you do not control: do not send it characters."""
    script("MethylEsterT3In", ENTERING, tunnel_obj="tun_cl_3").encode("ascii")


def test_it_says_what_is_missing_before_failing_thirty_times():
    """Run without the session and every line fails on its own, each describing a symptom."""
    text = script("AcidT3In", ENTERING, tunnel_obj="tun_cl_3")
    assert "cmd.get_object_list()" in text
    assert "session that is not loaded" in text
    guard = text[:text.index("delete snap_")]
    assert "AcidT3In" in guard and "tun_cl_3" in guard, "the check must name what it needs"


def test_the_guard_survives_having_no_tunnel_object():
    text = script("AcidT3In", ENTERING)
    assert "_missing" in text
    text.encode("ascii")
