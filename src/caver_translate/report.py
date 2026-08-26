"""What comes out: two tables and one page.

The tables are for working with -- open them in anything. The page is for reading: every profile
drawn, sorted so the easiest routes are at the top, with the caveats attached to the rows they
apply to rather than buried in a methods paragraph.

No plotting library. Profiles are drawn as inline SVG paths, so the report stays a single file that
opens anywhere and needs nothing installed.
"""
from __future__ import annotations

import csv
import html
from pathlib import Path

from .metrics import coverage, evaluate

FLAG_TEXT = {
    "short_tunnel": "tunnel shorter than 2 A: this is the mouth of the pocket, not a route through "
                    "the protein. It scores well because there is no distance to cross.",
    "positive_surface": "energy at the mouth is positive: the ligand already clashes there, so "
                        "dE_BS looks favourable only because a positive number was subtracted.",
    "direction_mismatch": "the tunnel widens the opposite way to what the file name claims. One of "
                          "the two is wrong; the radius was trusted.",
    "orientation_from_name": "constant radius, so which end is the binding site was taken from the "
                             "file name.",
    "lower_bound_only": "no upper-bound trajectory: the lower bound can pass through "
                        "discontinuities and understate a barrier.",
    "upper_bound_failed": "the upper bound was calculated and did not converge: with its rotation "
                          "constrained the ligand does not get past one of the discs. A result, "
                          "not a gap.",
    "failed": "no profile in the archive. CaverWeb keeps no log for a failed combination.",
}

COLUMNS = ["receptor", "ligand", "tunnel", "direction", "E_surface", "E_bound", "E_max",
           "Ea", "dE_BS", "n_discs", "span_A", "tunnel_length_A", "bottleneck_radius_A",
           "curvature", "priority", "flags", "source"]

TUNNEL_COLUMNS = ["receptor", "tunnel", "bottleneck_radius", "length", "curvature", "priority"]


def rows(tunnels, jobs) -> list:
    """One record per calculation, with the tunnel's geometry joined onto it."""
    by_id = {(t.receptor, t.tunnel): t for t in tunnels}
    out = []
    for job in jobs:
        tunnel = by_id.get((job.receptor, job.tunnel))
        m = evaluate(job, tunnel)
        out.append({
            "receptor": job.receptor, "ligand": job.ligand, "tunnel": job.tunnel,
            "direction": job.direction,
            "E_surface": m.energy_surface, "E_bound": m.energy_bound, "E_max": m.energy_max,
            "Ea": m.activation, "dE_BS": m.delta_bs,
            "n_discs": m.n_discs, "span_A": round(m.span, 2) if m.span is not None else None,
            "tunnel_length_A": tunnel.length if tunnel else None,
            "bottleneck_radius_A": tunnel.bottleneck_radius if tunnel else None,
            "curvature": tunnel.curvature if tunnel else None,
            "priority": tunnel.priority if tunnel else None,
            "flags": " ".join(m.flags), "source": job.source,
        })
    out.sort(key=lambda r: (r["receptor"], str(r["ligand"]), r["tunnel"] or 0, r["direction"] or ""))
    return out


def by_route(records) -> list:
    """One row per receptor, ligand and tunnel, with the two directions side by side.

    A route is the thing being chosen between, and it was calculated twice: entering and leaving.
    Read one direction at a time and the same tunnel appears in two places in the ranking, which
    is the one comparison nobody wants to make.

    E_surface and dE_BS are taken from the entering run when there is one. They describe the two
    ends of the route and both runs measure the same two ends, so the pair only differs by the
    noise of two separate dockings.
    """
    routes = {}
    for r in records:
        if r["tunnel"] is None:
            continue
        key = (r["receptor"], str(r["ligand"]), r["tunnel"])
        route = routes.setdefault(key, {
            "receptor": r["receptor"], "ligand": str(r["ligand"]), "tunnel": r["tunnel"],
            "ea_in": None, "ea_out": None, "dE_BS": None, "E_surface": None,
            "length": r["tunnel_length_A"], "neck": r["bottleneck_radius_A"], "flags": set(),
        })
        route["ea_in" if r["direction"] == "in" else "ea_out"] = r["Ea"]
        if r["direction"] == "in" or route["E_surface"] is None:
            route["dE_BS"], route["E_surface"] = r["dE_BS"], r["E_surface"]
        route["flags"].update(r["flags"].split() if r["flags"] else [])
    return list(routes.values())


def _easiest_first(routes) -> list:
    """Least resistance first, with the routes that were never calculated entering at the end."""
    return sorted(routes, key=lambda r: (r["ea_in"] is None, r["ea_in"] if r["ea_in"] is not None
                                         else (r["ea_out"] if r["ea_out"] is not None else 0.0)))


def tunnel_rows(tunnels) -> list:
    return [{"receptor": t.receptor, "tunnel": t.tunnel,
             "bottleneck_radius": t.bottleneck_radius, "length": t.length,
             "curvature": t.curvature, "priority": t.priority} for t in tunnels]


def write_csv(path, records, columns) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    return path


def spark(profile, width=250, height=52) -> str:
    """The energy profile as an SVG path, surface on the left, the barrier marked.

    Drawn by hand rather than with a plotting library so that the report has no dependency and
    stays one file: eighty of these in a folder of PNGs is a worse deliverable than one page.
    """
    if len(profile) < 2:
        return ""
    xs = [p.distance for p in profile]
    ys = [p.energy_lb for p in profile]
    dx = (max(xs) - min(xs)) or 1.0
    dy = (max(ys) - min(ys)) or 1.0
    pts = [(4 + (x - min(xs)) / dx * (width - 8), height - 6 - (y - min(ys)) / dy * (height - 14))
           for x, y in zip(xs, ys)]
    path_d = "M" + " L".join("%.1f,%.1f" % p for p in pts)
    peak = min(pts, key=lambda p: p[1])
    return ('<svg viewBox="0 0 %d %d" width="%d" height="%d">'
            '<path d="%s" fill="none" stroke="currentColor" stroke-width="1.6"/>'
            '<circle cx="%.1f" cy="%.1f" r="2.6" fill="#d33"></circle></svg>'
            % (width, height, width, height, path_d, peak[0], peak[1]))


STYLE = """
<style>
 body{font:14px/1.55 system-ui,-apple-system,sans-serif;margin:2rem auto;max-width:1150px;
      padding:0 1rem;color:#1b2a2e}
 h1{margin:0 0 .2rem} h2{margin:2rem 0 .3rem}
 .sub{color:#5a6b70;margin:0 0 1rem}
 table{border-collapse:collapse;width:100%;margin:.8rem 0}
 th,td{padding:.4rem .5rem;border-bottom:1px solid #e3e8ea;text-align:right;white-space:nowrap}
 th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}
 th{background:#f4f7f8;font-weight:600}
 .flag{color:#a4562a;font-size:12px;white-space:normal}
 .note{background:#fff8ee;border-left:3px solid #e0a15c;padding:.6rem .9rem;margin:.8rem 0;
       font-size:13px}
 svg{color:#2e7d8f}
 code{background:#f1f5f6;padding:.1rem .3rem;border-radius:3px}
</style>
"""


def _num(value) -> str:
    return "" if value is None else "%.2f" % value


def _route_table(routes, group, ranked_by, headers) -> str:
    """The routes grouped one way and ranked within each group, as html lines.

    The group's own columns are printed once and left blank on the rows below, so that what varies
    down the table is the ranking and not the two things being held fixed.
    """
    out = ["<table><tr><th>%s</th><th>%s</th><th>%s</th><th>Ea in</th><th>Ea out</th>"
           "<th>dE_BS</th><th>E_surface</th><th>Length A</th><th>Neck A</th><th>Notes</th></tr>"
           % tuple(html.escape(h) for h in headers)]
    keyed = {}
    for r in routes:
        keyed.setdefault(tuple(r[g] for g in group), []).append(r)
    for key in sorted(keyed, key=lambda k: tuple(str(part) for part in k)):
        for i, r in enumerate(_easiest_first(keyed[key])):
            head = ["", ""] if i else [html.escape(str(part)) for part in key]
            # dE_BS says nothing when the mouth is a clash, so it is not printed as though it did.
            bs = _num(r["dE_BS"])
            if "positive_surface" in r["flags"]:
                bs = '<span class="flag">%s (mouth clashes)</span>' % bs
            notes = ('<span class="flag">%s</span>' % html.escape(" / ".join(sorted(r["flags"])))
                     if r["flags"] else "")
            out.append("<tr><td>%s</td><td>%s</td><td>%s</td><td><strong>%s</strong></td>"
                       "<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
                       % (head[0], head[1], html.escape(str(r[ranked_by])),
                          _num(r["ea_in"]), _num(r["ea_out"]), bs, _num(r["E_surface"]),
                          _num(r["length"]), _num(r["neck"]), notes))
    out.append("</table>")
    return out


def write_html(path, tunnels, jobs) -> Path:
    records = rows(tunnels, jobs)
    by_source = {j.source: j for j in jobs}
    cov = coverage(jobs)
    ranked = sorted((r for r in records if r["Ea"] is not None), key=lambda r: r["Ea"])

    out = ['<!doctype html><meta charset="utf-8"><title>Ligand transport</title>', STYLE,
           "<h1>Ligand transport</h1>",
           '<p class="sub">%d calculations, %d tunnels, %d of %d combinations present</p>'
           % (len(records), len(tunnels), cov["present"], cov["expected"])]

    if cov["missing"] or cov["duplicated"]:
        out.append('<div class="note"><strong>Gaps.</strong> A CaverWeb combination that fails '
                   "leaves no log behind, so a missing row is the only trace that it was "
                   "attempted.<ul>")
        for receptor, ligand, tunnel, direction in cov["missing"][:40]:
            out.append("<li>missing: %s, %s, tunnel %s, %s</li>"
                       % (html.escape(receptor), html.escape(str(ligand)), tunnel, direction))
        for (receptor, ligand, tunnel, direction), sources in cov["duplicated"]:
            out.append("<li>two archives claim %s, %s, tunnel %s, %s: %s</li>"
                       % (html.escape(receptor), html.escape(str(ligand)), tunnel, direction,
                          html.escape(", ".join(sources))))
        out.append("</ul></div>")

    out.append("<h2>Easiest route first</h2>")
    out.append('<p class="sub">Ranked by <code>Ea</code>, what entering costs. Read it beside the '
               "tunnel length: a tunnel with nothing to cross always wins.</p>")
    out.append("<table><tr><th>Receptor</th><th>Ligand</th><th>Tunnel</th><th>Dir</th>"
               "<th>Ea</th><th>dE_BS</th><th>E_surface</th><th>E_max</th><th>E_bound</th>"
               "<th>Length A</th><th>Neck A</th><th>Profile</th></tr>")
    for r in ranked:
        job = by_source.get(r["source"])
        flags = r["flags"].split() if r["flags"] else []
        mark = '<div class="flag">%s</div>' % html.escape(" / ".join(flags)) if flags else ""
        out.append(
            "<tr><td>%s</td><td>%s%s</td><td>%s</td><td>%s</td>"
            "<td><strong>%s</strong></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
            "<td>%s</td><td>%s</td><td>%s</td></tr>"
            % (html.escape(r["receptor"]), html.escape(str(r["ligand"])), mark,
               r["tunnel"], r["direction"], _num(r["Ea"]), _num(r["dE_BS"]),
               _num(r["E_surface"]), _num(r["E_max"]), _num(r["E_bound"]),
               _num(r["tunnel_length_A"]), _num(r["bottleneck_radius_A"]),
               spark(job.profile) if job else ""))
    out.append("</table>")

    routes = by_route(records)

    out.append("<h2>Which tunnel each compound prefers</h2>")
    out.append('<p class="sub">Every route one compound has, least resistance first. '
               "<code>Ea</code> is what entering costs, so a low one is a way in and a high one is "
               "a real obstacle. <code>dE_BS</code> says how much better the site is than the "
               "surface, and it is only worth reading once <code>E_surface</code> is negative: a "
               "positive mouth means the ligand already clashes there, and subtracting a positive "
               "number makes the binding look excellent when nothing was measured. Read "
               "<code>Ea</code> beside the length as well; a tunnel with nothing to cross has "
               "nothing to cost.</p>")
    out += _route_table(routes, ("receptor", "ligand"), "tunnel",
                        ["Receptor", "Compound", "Tunnel"])

    out.append("<h2>Which compounds prefer each tunnel</h2>")
    out.append('<p class="sub">The same routes, read the other way: one tunnel and the compounds '
               "that get through it most easily.</p>")
    out += _route_table(routes, ("receptor", "tunnel"), "ligand",
                        ["Receptor", "Tunnel", "Compound"])

    used = sorted({f for r in records for f in (r["flags"].split() if r["flags"] else [])})
    if used:
        out.append("<h2>What the marks mean</h2><ul>")
        for flag in used:
            out.append("<li><strong>%s</strong> - %s</li>"
                       % (html.escape(flag), html.escape(FLAG_TEXT.get(flag, ""))))
        out.append("</ul>")

    out.append('<div class="note">These energies compare, they do not measure. CaverDock reports '
               "approximate docking energies along a path, not binding free energies, and the "
               "receptor is rigid unless flexible side chains were enabled.</div>")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out), encoding="utf-8")
    return path
