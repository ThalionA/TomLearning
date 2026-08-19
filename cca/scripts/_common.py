"""Shared plumbing for every script in ``cca/scripts/`` — import it FIRST.

    from _common import config, PAIRS, PAIR_NAMES, RES, results_path, fs_suffix, ...

What it does (each of these was copy-pasted into 20–80 scripts before 2026-08-19,
audit/REPORT_2026-08-19.md §3):

* puts ``cca/src`` on ``sys.path`` so ``tom_cca`` imports without an install
  (``pyproject.toml`` + ``uv pip install -e .`` makes this a no-op);
* re-exports ``config`` and the canonical locations ``RES`` / ``FIGS``;
* ``PAIRS`` (tuples) / ``PAIR_NAMES`` (strings) in the temporal-arm panel order;
* the FS-condition suffix in ONE idiom (``fs_suffix``, ``fs_from_argv``);
* ``results_path`` / ``results_stem`` — the ``<stem>_bin<ms><fs>.csv`` convention;
* ``temporal_parser`` — the common driver flags with the canonical defaults from
  ``config.TEMPORAL`` (``--bin-ms 10 --smooth-ms 2.5 --include-fs --max-samples 0 ...``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from tom_cca import config  # noqa: E402

RES = config.RESULTS_DIR
FIGS = config.FIGURES_DIR
PAIRS = list(config.PAIRS_TEMPORAL)       # [("CA1","RSC"), ...] panel order
PAIR_NAMES = list(config.PAIR_NAMES)      # ["CA1-RSC", ...]
TEMPORAL = config.TEMPORAL
EPOCHS = list(config.EPOCH_NAMES)         # ["naive", "intermediate", "expert"]
FS_SUFFIX = "_fsincl"
# The two co-primary FS conditions every analyze/figs script loops over: (suffix, label).
FS_CONDITIONS = [("", "FS-excluded"), (FS_SUFFIX, "FS-included")]


def fs_suffix(include_fs: bool) -> str:
    """Filename suffix for the FS condition: ``""`` (FS excluded, primary) or ``"_fsincl"``."""
    return FS_SUFFIX if include_fs else ""


def fs_from_argv(argv=None) -> str:
    """The reader idiom: ``python figs_x.py fsincl`` → ``"_fsincl"``, else ``""``.
    (Every analyze/figs script defaults to FS-EXCLUDED.)"""
    argv = sys.argv[1:] if argv is None else list(argv)
    return FS_SUFFIX if "fsincl" in argv else ""


def results_stem(stem: str, *, fs: str | bool = "", bin_ms: int | None = TEMPORAL.bin_ms) -> str:
    """``"<stem>_bin<bin_ms><fs>"`` — ``fs`` may be the suffix string or a bool."""
    suf = fs_suffix(fs) if isinstance(fs, bool) else fs
    b = "" if bin_ms is None else f"_bin{int(bin_ms)}"
    return f"{stem}{b}{suf}"


def results_path(stem: str, *, fs: str | bool = "", bin_ms: int | None = TEMPORAL.bin_ms,
                 ext: str = "csv") -> Path:
    """``RES / "<stem>_bin<bin_ms><fs>.<ext>"`` — the one place the convention lives."""
    return RES / f"{results_stem(stem, fs=fs, bin_ms=bin_ms)}.{ext}"


def fig_name(stem: str, *, fs: str | bool = "", bin_ms: int | None = TEMPORAL.bin_ms) -> str:
    """``"HCV1_<stem>_<fsexcl|fsincl>_bin<ms>"`` — the temporal figure stem convention."""
    suf = fs_suffix(fs) if isinstance(fs, bool) else fs
    fs_tag = "fsincl" if suf == FS_SUFFIX else "fsexcl"
    b = "" if bin_ms is None else f"_bin{int(bin_ms)}"
    return f"{config.FIG_PREFIX}_{stem}_{fs_tag}{b}"


def temporal_parser(description: str = "", *, cap: bool = True, shuffles: bool = True,
                    max_lag: bool = True, animals: bool = True) -> argparse.ArgumentParser:
    """Common driver flags with the canonical ``config.TEMPORAL`` defaults."""
    t = TEMPORAL
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--bin-ms", type=int, default=t.bin_ms, help=f"bin width ms (default {t.bin_ms})")
    p.add_argument("--smooth-ms", type=float, default=t.smooth_ms,
                   help=f"Gaussian s.d. (ms) on the 1 ms spike train before binning "
                        f"(default {t.smooth_ms}; 0 = raw counts)")
    p.add_argument("--include-fs", action="store_true",
                   help="include fast-spiking units (writes the _fsincl file; both conditions are co-primary)")
    p.add_argument("--out", default="", help="output CSV stem override")
    if max_lag:
        p.add_argument("--max-lag", type=int, default=t.max_lag_bins,
                       help=f"widest lag in bins; curve spans -max_lag..+max_lag (default {t.max_lag_bins})")
    if shuffles:
        p.add_argument("--n-shuffles", type=int, default=t.n_shuffles,
                       help=f"circular-shift surrogates (default {t.n_shuffles}; the empirical p floor "
                            f"is 1/(n+1) and BH over {t.fdr_dims} dims needs > {t.fdr_dims * 20 - 1})")
    if cap:
        p.add_argument("--max-samples", type=int, default=0,
                       help="cap on running bins per animal, consecutive within-trial blocks from the "
                            "EARLIEST trials (0 = no cap, the default; 12000 was the pre-2026-08-17 default)")
        p.add_argument("--max-per-trial", type=int, default=0,
                       help="max consecutive bins kept per trial when capping (0 = whole trial)")
    if animals:
        p.add_argument("--animals", default="",
                       help="comma-separated animal ids to run (smoke tests / timing)")
    return p


def animals_filter(arg: str):
    """``--animals "3,7"`` → {3, 7}; empty → None (all)."""
    return {int(x) for x in arg.split(",") if x.strip()} if arg else None


def cfg_from_args(args):
    """``config.DEFAULT`` with the driver's bin width / smoothing / FS flag applied."""
    import dataclasses
    return dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                               exclude_fast_spiking=not args.include_fs,
                               gaussian_sd_ms=args.smooth_ms)


def cap_desc(args) -> str:
    ms, mpt = getattr(args, "max_samples", 0), getattr(args, "max_per_trial", 0)
    if not ms and not mpt:
        return "NONE (all running bins)"
    return f"{ms or 'inf'} bins, <= {mpt or 'all'} per trial, EARLIEST trials first"
