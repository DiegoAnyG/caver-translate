"""Grouping the calculations into routes, so that a compound can be read one at a time.

A route was calculated twice, entering and leaving. Ranked as two separate rows the same tunnel
appears twice in its own ranking and is compared with itself, which answers nothing.
"""
from caver_translate.report import by_route, _easiest_first


def record(ligand, tunnel, direction, ea, dbs, surface, flags=""):
    return {"receptor": "R", "ligand": ligand, "tunnel": tunnel, "direction": direction,
            "Ea": ea, "dE_BS": dbs, "E_surface": surface, "tunnel_length_A": 5.0,
            "bottleneck_radius_A": 2.0, "flags": flags}


def test_the_two_directions_land_on_one_route():
    routes = by_route([record("benzo", 1, "in", 0.2, -1.0, -4.9),
                       record("benzo", 1, "out", 0.5, -1.1, -4.8)])
    assert len(routes) == 1
    assert (routes[0]["ea_in"], routes[0]["ea_out"]) == (0.2, 0.5)
    # Both runs measure the same two ends, and the entering one is the one the ranking is read
    # along, so its numbers are the ones shown.
    assert routes[0]["E_surface"] == -4.9


def test_a_route_calculated_one_way_only_keeps_what_it_has():
    routes = by_route([record("pent", 2, "out", 0.5, -1.0, -4.0)])
    assert routes[0]["ea_in"] is None and routes[0]["ea_out"] == 0.5
    assert routes[0]["E_surface"] == -4.0


def test_flags_from_either_direction_stay_on_the_route():
    routes = by_route([record("et", 4, "in", 3.7, -8.2, 2.7, "positive_surface"),
                       record("et", 4, "out", 3.4, -8.4, 2.9, "short_tunnel")])
    assert routes[0]["flags"] == {"positive_surface", "short_tunnel"}


def test_a_failed_archive_has_no_route():
    assert by_route([record("x", None, None, None, None, None, "failed")]) == []


def test_least_resistance_first_and_the_uncalculated_last():
    routes = by_route([record("a", 1, "in", 1.5, -1.0, -4.0),
                       record("b", 2, "in", 0.2, -1.0, -4.0),
                       record("c", 3, "out", 0.1, -1.0, -4.0)])
    assert [r["tunnel"] for r in _easiest_first(routes)] == [2, 1, 3]
