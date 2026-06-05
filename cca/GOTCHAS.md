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
