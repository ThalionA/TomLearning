"""Stage 3 driver for the Tom sweep.

Communication-subspace membership and principal-angle reorientation for every
config in the spatial sweep. No surrogate loop, so it runs fast (a few seconds
per config); resumable -- configs whose results/stage3_<tag>.pkl already exists
are skipped unless ``--fresh``.

Run:  python scripts/run_stage3.py [--fresh] [--data-dir PATH]
"""

from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tom_cca import config, dataio, pipeline, stage3, sweep  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sweep", choices=("spatial",), default="spatial")
    p.add_argument("--fresh", action="store_true",
                   help="recompute configs even if a result file exists")
    p.add_argument("--data-dir", type=str, default=None,
                   help="override config.DATA_DIR (where TF*_export.mat lives)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    grid = sweep.build_sweep(args.sweep)
    print(f"Sweep '{args.sweep}': Stage 3 for {len(grid)} configs.\n")
    animals = dataio.load_animals(args.data_dir,
                                  spatial_field=config.DEFAULT.spatial_field)
    behaviour = dataio._read_behaviour_file(
        (Path(args.data_dir) if args.data_dir else config.DATA_DIR)
        / config.LEARNING_FILE
    )
    for tag, cfg in grid:
        out = config.RESULTS_DIR / f"stage3_{tag}.pkl"
        if out.exists() and not args.fresh:
            print(f"  skip '{tag}' (done)")
            continue
        run_variant(animals, behaviour, cfg, tag)
    print("\nStage 3 sweep done.")


def run_variant(animals, behaviour, cfg, tag: str) -> None:
    entries = dataio.classify_cohort(animals, cfg, behaviour)
    results: list[stage3.PairSubspace] = []
    t0 = time.time()
    for animal in animals:
        if animal.animal_id not in entries:        # learners only
            continue
        entry = entries[animal.animal_id]
        for area_x, area_y in config.PAIRS:
            prepared = pipeline.prepare_pair(animal, area_x, area_y, entry, cfg)
            if isinstance(prepared, pipeline.PreparedPair):
                results.append(stage3.analyse_subspace(prepared, cfg))
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = config.RESULTS_DIR / f"stage3_{tag}.pkl"
    with open(out, "wb") as fh:
        pickle.dump({"results": results, "cfg": cfg}, fh)
    print(f"  [{tag}] ({pipeline.config_label(cfg)}) {len(results)} pairs "
          f"in {time.time() - t0:.1f}s -> {out.name}")


if __name__ == "__main__":
    main()
