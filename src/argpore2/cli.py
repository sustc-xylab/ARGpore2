"""Command-line entry point. The analysis itself remains argpore.sh."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from argpore2 import __version__

FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".fna", ".ffn"}


def repo_root() -> Path:
    env = os.environ.get("ARGPORE_HOME", "").strip()
    if env:
        root = Path(env).expanduser().resolve()
        if (root / "argpore.sh").is_file():
            return root
        raise FileNotFoundError(f"ARGPORE_HOME has no argpore.sh: {root}")
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "argpore.sh").is_file() and (parent / "check-setup.sh").is_file():
            return parent
    raise FileNotFoundError(
        "ARGpore2 repository not found. Install from the clone with "
        "'pip install -e .' or set ARGPORE_HOME."
    )


def conda_env_bin(env_name: str) -> Path | None:
    prefix = os.environ.get("CONDA_PREFIX", "")
    candidates: list[Path] = []
    if prefix:
        pp = Path(prefix)
        if pp.parent.name == "envs":
            candidates.append(pp.parent.parent)
        else:
            candidates.append(pp)
    conda = shutil_which("conda")
    if conda:
        exe = Path(conda).resolve()
        if exe.parent.name in {"bin", "condabin"}:
            candidates.append(exe.parent.parent)
    for root in candidates:
        bindir = root / "envs" / env_name / "bin"
        if bindir.is_dir():
            return bindir
    return None


def shutil_which(name: str) -> str | None:
    from shutil import which

    return which(name)


def runtime_env(root: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PERL5LIB", None)
    env["ARGPORE_HOME"] = str(root)
    r_env = os.environ.get("ARGPORE_R_ENV", "argpore2-r")
    bindir = conda_env_bin(r_env)
    if bindir and (bindir / "Rscript").is_file():
        env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
    return env


def run_script(root: Path, script: str, args: list[str]) -> int:
    path = root / script
    if not path.is_file():
        print(f"[ERROR] missing {path}", file=sys.stderr)
        return 1
    completed = subprocess.run(
        ["bash", str(path), *args],
        cwd=root,
        env=runtime_env(root),
    )
    return completed.returncode


def cmd_check() -> int:
    return run_script(repo_root(), "check-setup.sh", [])


def cmd_setup_tools() -> int:
    return run_script(repo_root(), "setup-tools.sh", [])


def cmd_setup_db() -> int:
    return run_script(repo_root(), "setup-db.sh", [])


def cmd_run(args: argparse.Namespace) -> int:
    fasta = Path(args.fasta).expanduser().resolve()
    if not fasta.is_file():
        print(f"[ERROR] input FASTA not found: {fasta}", file=sys.stderr)
        return 1
    if fasta.suffix.lower() not in FASTA_SUFFIXES:
        accepted = " ".join(sorted(FASTA_SUFFIXES))
        print(
            f"[ERROR] unsupported suffix {fasta.suffix!r}. Accepted: {accepted}",
            file=sys.stderr,
        )
        return 1
    root = repo_root()
    cmd = [
        "bash",
        str(root / "argpore.sh"),
        "-f",
        str(fasta),
        "-t",
        str(args.threads),
        "-o",
        str(Path(args.outdir).expanduser().resolve()),
        "-l",
        str(args.length),
        "-s",
        str(args.similarity),
    ]
    completed = subprocess.run(cmd, env=runtime_env(root))
    return completed.returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="argpore2",
        description=(
            "Identify ARGs and their hosts from nanopore reads. "
            "Maintenance commands: argpore2 check, setup-tools, setup-db."
        ),
    )
    parser.add_argument("--version", action="version", version=f"argpore2 {__version__}")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA (.fa, .fasta, .fas, .fna, .ffn)")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-o", "--outdir", required=True, help="Directory for all outputs")
    parser.add_argument("-l", "--length", default="0.9", help="Alignment length cutoff")
    parser.add_argument("-s", "--similarity", default="0.7", help="Similarity cutoff")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["check"]:
        return cmd_check()
    if argv[:1] == ["setup-tools"]:
        return cmd_setup_tools()
    if argv[:1] == ["setup-db"]:
        return cmd_setup_db()
    return cmd_run(build_parser().parse_args(argv))
