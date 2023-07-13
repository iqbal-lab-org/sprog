import os
import pytest

from sprog import utils


def test_load_and_write_json():
    filename = "tmp.test.json"
    utils.syscall(f"rm -f {filename}")
    data = {"foo": "bar"}
    utils.write_json(data, filename)
    got = utils.load_json(filename)
    assert got == data
    os.unlink(filename)
