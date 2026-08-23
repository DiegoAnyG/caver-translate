"""Read CaverWeb downloads without opening eighty archives by hand.

CaverWeb runs one job per ligand, tunnel and direction, and hands each one back as a zip named
after a hash. A modest study is eighty of them. This turns that folder into a table of the numbers
that get reported -- what entering costs, where the barrier is, how good the destination is -- with
the things that mislead those numbers marked on the rows they affect.

    from caver_translate import scan, rows
    tunnels, jobs = scan("CaverWEB/")
    table = rows(tunnels, jobs)
"""
from .metrics import Metrics, coverage, evaluate      # noqa: F401
from .parse import Job, Point, Tunnel, parse_job, parse_tunnels, scan   # noqa: F401
from .report import rows, tunnel_rows, write_csv, write_html            # noqa: F401

try:
    from importlib.metadata import version as _version

    __version__ = _version("caver-translate")
except Exception:                                     # a source tree with nothing installed
    __version__ = "0.1.0"

__all__ = ["scan", "rows", "evaluate", "coverage", "write_csv", "write_html",
           "Job", "Point", "Tunnel", "Metrics", "parse_job", "parse_tunnels", "tunnel_rows"]
