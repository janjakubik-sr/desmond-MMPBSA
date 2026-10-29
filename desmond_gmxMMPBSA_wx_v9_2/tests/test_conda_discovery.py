#!/usr/bin/env python3
from __future__ import annotations

import ast
import os
import re
import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "desmond_gmxmmpbsa_wx.py"
WANTED = {
    "_resolve_executable_candidate",
    "probe_conda",
    "discover_conda",
    "resolve_conda_executable",
}


def load_functions():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in WANTED
    ]
    namespace = {
        "os": os,
        "re": re,
        "shutil": shutil,
        "subprocess": subprocess,
        "Path": Path,
        "lru_cache": lru_cache,
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), namespace)
    return namespace


def main() -> None:
    functions = load_functions()
    with tempfile.TemporaryDirectory() as directory:
        home = Path(directory)
        bad = home / "bin" / "conda"
        good = home / "miniforge3" / "bin" / "conda"
        bad.parent.mkdir(parents=True)
        good.parent.mkdir(parents=True)

        bad.write_text(
            '#!/bin/sh\necho "ArgumentError: activate does not accept more than one argument" >&2\nexit 2\n'
        )
        good.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "--version" ]; then echo "conda 25.5.1"; exit 0; fi\n'
            'if [ "$1" = "run" ] && [ "$2" = "--help" ]; then echo "run help"; exit 0; fi\n'
            'exit 2\n'
        )
        bad.chmod(0o755)
        good.chmod(0o755)

        old = {key: os.environ.get(key) for key in ("HOME", "PATH", "CONDA_EXE")}
        try:
            os.environ["HOME"] = str(home)
            os.environ["PATH"] = str(bad.parent) + os.pathsep + (old["PATH"] or "")
            os.environ["CONDA_EXE"] = str(bad)
            functions["probe_conda"].cache_clear()
            functions["discover_conda"].cache_clear()
            assert functions["discover_conda"]() == str(good.resolve())
            assert functions["resolve_conda_executable"](str(bad)) == str(good.resolve())
        finally:
            for key, value in old.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value

    print("Conda discovery test passed")


if __name__ == "__main__":
    main()
