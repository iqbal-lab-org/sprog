from pkg_resources import get_distribution

try:
    __version__ = get_distribution("sprog").version
except:
    __version__ = "local"


__all__ = [
    "bifrost",
    "fastani",
    "genome_download",
    "gtdb",
    "mash",
    "sample",
    "sample_set",
    "tasks",
    "tree",
    "utils",
]

from sprog import *
