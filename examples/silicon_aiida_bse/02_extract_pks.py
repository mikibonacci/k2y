#!/usr/bin/env python
"""
02_extract_pks.py – extract AiiDA PKs from a Koopmans step_data.pkl
=====================================================================

After `koopmans run --engine aiida` finishes (step A), the AiiDA engine writes
a step_data.pkl that records the AiiDA PK of every submitted workchain, keyed
by the step's unique ID.

This script reads that file and prints the two PKs required by the next steps:

  KOOPMANS_NSCF_PK  – the PwBaseWorkChain for the **final** (dense-grid) NSCF.
                      This node carries the structure, k-mesh and pw.x save that
                      was converted to a Yambo SAVE directory.  Feed this into
                      01_p2y.py.

  KOOPMANS_HAM_PK   – the KcwCalculation (ham) that computed the pKI Koopmans
                      eigenvalues.  Feed this into 03_bse_submit.py.

Usage
-----
  python 02_extract_pks.py [path/to/step_data.pkl]

  If no path is given, looks for step_data.pkl in the current directory.
"""

import sys
import pathlib

# ---------------------------------------------------------------------------
# Locate the pkl file
# ---------------------------------------------------------------------------
if len(sys.argv) > 1:
    pkl_path = pathlib.Path(sys.argv[1])
else:
    pkl_path = pathlib.Path("step_data.pkl")

if not pkl_path.exists():
    print(f"ERROR: '{pkl_path}' not found.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
try:
    import dill as pickle
except ImportError:
    import pickle

with open(pkl_path, "rb") as fh:
    data = pickle.load(fh)

steps = data.get("steps", {})

# ---------------------------------------------------------------------------
# Print full step listing
# ---------------------------------------------------------------------------
print("=" * 68)
print(f"  Full step_data: {pkl_path}")
print("=" * 68)
for k, v in steps.items():
    print(f"  {k}: {v}")
print()

# ---------------------------------------------------------------------------
# Identify key steps
# ---------------------------------------------------------------------------
# Step naming convention produced by the Koopmans ASE workflow (method=dfpt):
#
#   01-koopmans-dfpt/
#     01-koopmans-dfpt-coarse/
#       01-wannierize/
#         02-nscf                        <- coarse NSCF  (NOT needed here)
#       02-kcw_wannier
#       03-compute-screening-via-dfpt/
#         01-kcw_screen
#     02-wannierize/
#       02-nscf                          <- FINAL dense-grid NSCF  ← NSCF_PK
#     03-kcw_wannier
#     04-kcw_ham                         ← HAM_PK
#
# We search for:
#   * the LAST step ending with "02-nscf" that does NOT belong to the coarse block
#   * the step ending with "kcw_ham"

nscf_candidates = {}
ham_pk = None

for uid, info in steps.items():
    wc = info.get("workchain")
    if wc is None:
        continue
    if uid.endswith("02-nscf"):
        nscf_candidates[uid] = wc
    if uid.endswith("kcw_ham"):
        ham_pk = wc

# Dense NSCF: in the "02-wannierize" block, not in "coarse"
dense_nscf_pk  = None
dense_nscf_uid = None
for uid in sorted(nscf_candidates):
    if "coarse" not in uid.lower():
        dense_nscf_pk  = nscf_candidates[uid]
        dense_nscf_uid = uid

# Fallback: take the last one
if dense_nscf_pk is None and nscf_candidates:
    dense_nscf_uid, dense_nscf_pk = sorted(nscf_candidates.items())[-1]

# ---------------------------------------------------------------------------
# Print results
# ---------------------------------------------------------------------------
print("=" * 68)
print("  PKs to use in 01_p2y.py and 03_bse_submit.py")
print("=" * 68)

if dense_nscf_pk is not None:
    print(f"  KOOPMANS_NSCF_PK = {dense_nscf_pk}   (step: {dense_nscf_uid})")
else:
    print("  KOOPMANS_NSCF_PK = NOT FOUND – check step names above")

if ham_pk is not None:
    print(f"  KOOPMANS_HAM_PK  = {ham_pk}")
else:
    print("  KOOPMANS_HAM_PK  = NOT FOUND – check step names above")

print()
print("Next steps:")
print("  1. Set KOOPMANS_NSCF_PK in 01_p2y.py and run it.")
print("  2. After 01_p2y.py finishes, set P2Y_WF_PK and KOOPMANS_HAM_PK in")
print("     03_bse_submit.py and run it.")
