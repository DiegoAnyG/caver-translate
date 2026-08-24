"""Choosing which trajectories to draw, without inventing a selection language for it.

The archive names carry the compound, the tunnel and the direction, so a shell glob already picks
any set a person could want. All this side has to do is read the tunnel and the hash back out of
the name, and expand a folder without dropping anything.
"""
from caver_translate.cli_pymol import TAIL, collect


def tail(name):
    m = TAIL.search(name)
    return None if m is None else (m["tunnel"], m["hash"])


def test_tunnel_and_hash_come_from_the_name():
    assert tail("benzo3inuh6v4fpkxc_results.zip") == ("3", "uh6v4fpkxc")
    assert tail("et1outihuhev2gt6_results.zip") == ("1", "ihuhev2gt6")
    assert tail("met3in4ywxjawqzf_results.zip") == ("3", "4ywxjawqzf")


def test_a_second_download_is_still_read():
    # "..._results (1).zip" is what a browser calls the second copy, and it is a whole
    # calculation: a name pattern anchored on the end of the file name would lose it.
    assert tail("benzo3inuh6v4fpkxc_results (1).zip") == ("3", "uh6v4fpkxc")


def test_an_unrecognised_name_is_not_guessed_at():
    assert tail("something_else.zip") is None


def test_a_folder_expands_to_its_archives(tmp_path):
    for name in ["b3inaaa_results.zip", "b3inaaa_results (1).zip", "notes.txt"]:
        (tmp_path / name).touch()
    assert [p.name for p in collect([tmp_path])] == ["b3inaaa_results (1).zip",
                                                     "b3inaaa_results.zip"]


def test_a_named_file_is_taken_as_given(tmp_path):
    f = tmp_path / "b3inaaa_results.zip"
    f.touch()
    assert collect([f]) == [f]


def test_a_wsl_path_is_printed_the_way_pymol_can_open_it():
    # PyMOL runs on Windows while this runs in WSL: /mnt/c/... after an @ is a file it cannot find.
    from pathlib import PurePosixPath
    from caver_translate.cli_pymol import for_pymol
    assert for_pymol(PurePosixPath("/mnt/c/Docs/poses/a.pml")) == "C:/Docs/poses/a.pml"
    assert for_pymol(PurePosixPath("/home/me/poses/a.pml")) == "/home/me/poses/a.pml"
