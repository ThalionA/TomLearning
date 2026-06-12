# GOTCHAS

One-line entries for non-obvious bugs, so they are not reintroduced.

- **CSV files use CRLF line endings.** The `csv` module writes `\r\n` by default even
  with `open(..., newline="")`. So `learning_changes_*.csv`, `landmark_prune_*.csv`, and
  `sweep_landmark_summary.csv` have a trailing `\r` on every field-10 value. A naive
  `awk -F, '$10=="True"'` matches **nothing** (it compares against `"True\r"`). Strip it
  first: `awk -F, '{gsub(/\r/,"")} ...'`, or read with pandas (handles it). Fixed going
  forward by passing `lineterminator="\n"` to the CSV writers (2026-06-05); already-written
  CSVs keep their CRLF until regenerated.

- **MATLAB `string`-type fields are unreadable by h5py.** `TF*_export.mat` region labels /
  learning points stored as MATLAB `string` cannot be read directly. Worked around with the
  one-off `scripts/export_cca_labels.m`, which writes `cca_labels.json` (schema
  `tom_cca_labels_v1`) that `dataio.py` reads as a companion file.

- **Judge landmark-config overfitting by `frac_cc_ge_099_*`, not `max_cc`.** A single
  saturated canonical dim pushes `max_cc` to 0.999 in otherwise-healthy configs. See
  `STATE.md` §4 and `src/tom_cca/prune_table.py`.

- **A single trial cannot re-fit pCCA/KCCA.** One trial ≈ 580 running 50 ms bins; a
  30-component fit needs ~50×30 ≈ 1500 samples. For per-trial resolution, PROJECT the
  trial onto a subspace fit on many other trials (`early_trials.reference_fit`), don't
  re-fit; fit-only metrics (weight-Gini, angles, #sig, KCCA) need ≥5-trial blocks. See
  `early_trials.py` and report §2.10c.

- **Don't nest `nohup … &` inside a `run_in_background` Bash call.** It double-backgrounds:
  the harness marks the *wrapper* complete (after the 2 s echo/sleep) while the real Python
  run keeps going untracked, so you never get a true completion signal. Launch the bare
  `python …` command with `run_in_background` (no `&`/nohup), or poll with `pgrep`.

- **The ~80° top-3 split-half "noise floor" is NOT a bug — it means the subspace is ~1-D.**
  `subspace_window.split_half_x` is the *max* principal angle over the top-3 canonical dims.
  When only CC1 carries a stable direction (dims 2–3 have cc≈0, i.e. no real structure), the
  max angle is dominated by the random noise dims → near-orthogonal. Verified: 3 genuinely
  shared dims → top-3 floor ~9°; 1 shared dim → CC1 floor ~7° but top-3 floor ~85°
  (`test_subspace_window.py::test_split_half_floor_tracks_true_dimensionality`). **Use the
  CC1-only floor (`split_half_x_cc1`) to judge dominant-direction stability**, never the
  top-3 max angle when n_sig is low. (report §3.4)

- **Running `run_trajectory` truncates its window CSV mid-run (`w` mode).** The driver opens
  `trajectory_w15_bin{25,50}{,_fsincl}.csv` in `"w"` and re-writes after each animal, so while a stage
  is in flight the file holds only the animals done so far. Regenerating trajectory figures against it
  mid-run silently builds them from 1–2 animals (caught 2026-06-12: the 23:56 fsexcl figures). **Only
  run `figs_report.py <traj-stem>` when the window CSV shows 16 animals**; to plot full data while a
  re-run is in flight, restore from git (`git show HEAD:cca/results/<stem>.csv > <stem>_safe.csv`) and
  point figs at the protected stem. The per-dim `_dims.csv` truncates the same way.

- **Too many `run_in_background` tasks → the harness evicts (kills) one.** Launching two batch
  runs as background tasks *plus* a watcher (3 tracked tasks) got a batch killed ~14 min in
  (2026-06-12) — looked like a crash but was eviction (RAM was 87 % free, not OOM). For long
  unattended compute, launch it **fully detached** so it isn't a harness task:
  `(nohup bash run.sh >/dev/null 2>&1 </dev/null &)` → reparents to launchd (PPID 1, new session),
  survives. Keep at most ~1 `run_in_background` watcher alongside.
