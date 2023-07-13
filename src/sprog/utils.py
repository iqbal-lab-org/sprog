import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys


def syscall(command, cwd=None, quiet=False):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        cwd=cwd,
    )
    logging.info(f"Return code: {completed_process.returncode}; {command}")
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise Exception("Error in system call. Cannot continue")

    if not quiet:
        logging.info(f"stdout:\n{completed_process.stdout.rstrip()}")
        logging.info(f"stderr:\n{completed_process.stderr.rstrip()}")
    return completed_process


def load_json(filename):
    with open(filename) as f:
        return json.load(f)


def write_json(data, filename):
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


def md5(filename):
    # see https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def copy_or_symlink(source, dest, symlink=False):
    if symlink:
        os.symlink(os.path.abspath(source), os.path.abspath(dest))
    else:
        shutil.copy(source, dest)
        if md5(source) != md5(dest):
            raise Exception(f"Error copying file (md5 mismatch): {source} -> {dest}")
