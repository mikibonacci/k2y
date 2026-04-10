# BSE@Koopmans with AiiDA – Silicon example

Full AiiDA-managed pipeline from Koopmans pKI-DFPT to a BSE optical spectrum
for bulk silicon (LDA), using k2y to map Koopmans eigenvalues into Yambo's QP
database (`ndb.QP`).

All calculations are submitted to a scheduler/HPC cluster via the AiiDA
framework.  Provenance (inputs, outputs, parameters) is recorded automatically
in the AiiDA database.

## Workflow overview

```
Step A  00_koopmans_aiida.sh   (koopmans run --engine aiida)
            pw.x SCF → pw.x NSCF → Wannier90 → kcw.x screening → kcw.x ham
            ↓ produces step_data.pkl

Step B  01_p2y.py              (aiida-yambo YamboWorkflow, INITIALISE=True)
            pw.x SCF → pw.x NSCF → p2y → yambo init
            ↓ produces Yambo SAVE directory on remote

Step C  02_extract_pks.py
            reads step_data.pkl, prints KOOPMANS_NSCF_PK + KOOPMANS_HAM_PK

Step D  03_bse_submit.py       (k2y + aiida-yambo YamboWorkflow)
            k2y.aiida.generate_kcw_qp_database → ndb.QP
            → YamboWorkflow BSE (reuses SAVE from step B)
```

## Files

| File | Description |
|---|---|
| `00_koopmans_aiida.sh` | Run Koopmans pKI via AiiDA engine; calls `02_extract_pks.py` at the end |
| `01_p2y.py` | Submit p2y YamboWorkflow (SAVE generation) |
| `02_extract_pks.py` | Extract `KOOPMANS_NSCF_PK` / `KOOPMANS_HAM_PK` from `step_data.pkl` |
| `03_bse_submit.py` | Generate `ndb.QP` via k2y (AiiDA calcfunction), submit BSE YamboWorkflow |
| `inputs/si_ki.json` | Koopmans pKI-DFPT input for bulk silicon |
| `inputs/aiida_engine.json` | AiiDA code labels and scheduler resources template |

## Prerequisites

- `koopmans`, `aiida-koopmans`, `aiida-yambo`, `k2y` installed
- AiiDA daemon running (`verdi daemon start`)
- Codes registered in AiiDA (labels must match `inputs/aiida_engine.json`)
- Pseudopotential family installed: `PseudoDojo/0.4/LDA/SR/standard/upf`

## Quick start

```bash
# 1. Edit inputs/aiida_engine.json to match your registered AiiDA code labels.
#    Edit 01_p2y.py and 03_bse_submit.py: set PW_CODE, PREPROCESSING_CODE, YAMBO_CODE.

# 2. Run Koopmans via AiiDA (blocking – submits workchains to the scheduler)
bash 00_koopmans_aiida.sh

# 3. Copy KOOPMANS_NSCF_PK from the output into 01_p2y.py, then run it
python 01_p2y.py

# 4. Once 01_p2y.py finishes, copy KOOPMANS_HAM_PK and P2Y_WF_PK into
#    03_bse_submit.py, then run it
python 03_bse_submit.py

# 5. Monitor progress
verdi process list -a
```

The BSE spectrum is available in the outputs of the submitted YamboWorkflow node.

## Documentation

A full step-by-step tutorial is available in the k2y documentation:

  docs/tutorials/silicon_bse_aiida.md
