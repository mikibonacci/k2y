#!/usr/bin/env bash
# =============================================================================
# Run a Koopmans pKI-DFPT workflow via the AiiDA engine – bulk silicon
# =============================================================================
#
# This script is step A of the AiiDA-based BSE@Koopmans workflow:
#
#   Step A  00_koopmans_aiida.sh   – run Koopmans (pw + W90 + kcw) via AiiDA
#   Step B  01_p2y.py              – submit p2y YamboWorkflow (SAVE generation)
#   Step C  02_extract_pks.py      – extract AiiDA PKs from step_data.pkl
#   Step D  03_bse_submit.py       – generate ndb.QP via k2y, submit BSE
#
# This script:
#   1. Copies si_ki.json and aiida_engine.json into a fresh working directory.
#   2. Runs `koopmans run` with --engine aiida and --engine_config.
#   3. After the run finishes it reads the generated step_data.pkl and prints
#      the PKs required by 01_p2y.py and 03_bse_submit.py.
#
# Usage
# -----
#   bash 00_koopmans_aiida.sh               # use defaults
#   WORKDIR=/some/other/path bash 00_koopmans_aiida.sh
#
# Requirements
# ------------
#   - koopmans installed and in PATH
#   - AiiDA daemon running       (verdi daemon start)
#   - Codes registered in AiiDA  (labels must match inputs/aiida_engine.json)
#   - Pseudopotentials installed (PseudoDojo/0.4/LDA/SR/standard/upf)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUTS="$SCRIPT_DIR/inputs"

# Working directory for the Koopmans run (created if it does not exist).
WORKDIR="${WORKDIR:-$SCRIPT_DIR/KOOPMANS_AIIDA}"

echo "================================================================"
echo "  Koopmans pKI (AiiDA engine) – Silicon"
echo "  Working dir: $WORKDIR"
echo "================================================================"
echo ""

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
for f in si_ki.json aiida_engine.json; do
    [[ -f "$INPUTS/$f" ]] || {
        echo "ERROR: $INPUTS/$f not found."
        echo "  Edit inputs/aiida_engine.json to match your registered AiiDA codes."
        exit 1
    }
done

# ---------------------------------------------------------------------------
# Prepare working directory
# ---------------------------------------------------------------------------
mkdir -p "$WORKDIR"
cp "$INPUTS/si_ki.json"        "$WORKDIR/"
cp "$INPUTS/aiida_engine.json" "$WORKDIR/"

cd "$WORKDIR"

# ---------------------------------------------------------------------------
# Check AiiDA daemon
# ---------------------------------------------------------------------------
if ! verdi daemon status 2>&1 | grep -q "running"; then
    echo "WARNING: AiiDA daemon does not appear to be running."
    echo "  Start it with:  verdi daemon start"
    echo "  Continuing anyway – processes will be queued."
fi

# ---------------------------------------------------------------------------
# Run Koopmans via AiiDA engine
# ---------------------------------------------------------------------------
# Flags:
#   -t  : print traceback on error
#   -l  : enable koopmans logging (koopmans.log)
#   --engine aiida         : use the AiiDA engine (submits to scheduler/cluster)
#   --engine_config ...    : JSON with code labels and scheduler resources
#
# The run is blocking: it polls until all AiiDA workchains finish.
# To run in background:  nohup bash 00_koopmans_aiida.sh &  tail -f nohup.out
# ---------------------------------------------------------------------------
echo "Running: koopmans run -t -l --engine aiida --engine_config aiida_engine.json si_ki.json"
koopmans run -t -l \
    --engine aiida \
    --engine_config aiida_engine.json \
    si_ki.json

echo ""
echo "--- Koopmans run finished ---"
echo ""

# ---------------------------------------------------------------------------
# Summarise the step_data.pkl
# ---------------------------------------------------------------------------
PKL="$WORKDIR/step_data.pkl"
if [[ -f "$PKL" ]]; then
    echo "=== step_data.pkl summary ==="
    aiida-koopmans explore "$PKL" 2>/dev/null || true
    echo ""
    echo "=== Extracting PKs for 01_p2y.py / 03_bse_submit.py ==="
    python "$SCRIPT_DIR/02_extract_pks.py" "$PKL"
else
    echo "WARNING: step_data.pkl not found at $PKL"
    echo "  The workflow may not have produced it yet."
    echo "  Run 02_extract_pks.py manually once it exists."
fi
