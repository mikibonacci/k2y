# Tutorial: Silicon BSE with Koopmans eigenvalues via AiiDA

This tutorial shows how to run the BSE@Koopmans workflow for bulk silicon
using [AiiDA](https://www.aiida.net) to manage all calculations.  Every step
is submitted to a scheduler/HPC cluster and full provenance is recorded in
the AiiDA database.

For the equivalent standalone (no-AiiDA) workflow see:
[Silicon BSE with Koopmans](silicon_bse.md).

The complete input files and scripts are in `examples/silicon_aiida_bse/`.

## What you need

Everything you need can be installed via `pip install .[aiida]`.

| Software | Purpose |
|---|---|
| `koopmans` + `aiida-koopmans` | Run the pKI-DFPT workflow via AiiDA |
| `aiida-yambo` | Submit Yambo calculations (p2y + BSE) via AiiDA |
| `k2y` | Map Koopmans eigenvalues → Yambo `ndb.QP` (with AiiDA provenance) |
| AiiDA daemon | Running and monitoring submitted workchains |

Pseudopotentials: `PseudoDojo/0.4/LDA/SR/standard/upf` (install with
`aiida-pseudo install pseudo-dojo`).

---

## Workflow overview

```
Step A  00_koopmans_aiida.sh
         koopmans run --engine aiida
         pw.x SCF → NSCF → Wannier90 → kcw.x screen → kcw.x ham
         └─ produces: KOOPMANS_AIIDA/step_data.pkl

Step B  01_p2y.py
         YamboWorkflow (INITIALISE=True)
         pw.x SCF → NSCF → p2y → yambo init
         └─ produces: Yambo SAVE on remote machine

Step C  02_extract_pks.py
         reads step_data.pkl
         └─ prints: KOOPMANS_NSCF_PK, KOOPMANS_HAM_PK

Step D  03_bse_submit.py
         k2y.aiida.generate_kcw_qp_database  (calcfunction)
         └─ ndb.QP stored as SinglefileData in AiiDA database
         YamboWorkflow (BSE, reuses SAVE from step B)
         └─ output: optical spectrum ε₂(ω)
```

---

## Before you start

### 1. Register AiiDA codes

Make sure `pw.x`, `p2y`, `yambo`, `projwfc.x`, `pw2wannier90.x`, `wannier90`,
and `kcw.x` are registered in AiiDA (see `verdi code list`).

### 2. Edit the engine configuration

Open `inputs/aiida_engine.json` and replace the placeholder labels with the
names of your registered codes:

```json
{
    "pw_code":           "pw@mycomputer",
    "projwfc_code":      "projwfc@mycomputer",
    "pw2wannier90_code": "pw2wannier90@mycomputer",
    "wannier90_code":    "wannier90@mycomputer",
    "kcw_code":          "kcw@mycomputer",
    ...
}
```

Adjust `num_mpiprocs_per_machine` and `max_wallclock_seconds` to match your
cluster.

### 3. Edit the Python scripts

At the top of `01_p2y.py` and `03_bse_submit.py`, set the code labels and
pseudo-family to match your AiiDA setup:

```python
PW_CODE            = "pw@mycomputer"
PREPROCESSING_CODE = "p2y@mycomputer"
YAMBO_CODE         = "yambo@mycomputer"
PSEUDO_FAMILY      = "PseudoDojo/0.4/LDA/SR/standard/upf"
```

### 4. Start the AiiDA daemon

```bash
verdi daemon start 4   # use as many workers as needed
```

---

## Step A — Koopmans pKI via AiiDA

```bash
bash 00_koopmans_aiida.sh
```

This runs the full Koopmans pKI-DFPT workflow:

1. `pw.x` SCF on a coarse 2×2×2 k-grid
2. `pw.x` NSCF on the same grid (nosym, extra bands) for Wannierisation
3. `wannier90` + `pw2wannier90` to build MLWFs
4. `kcw.x` screening (coarse grid → renormalise)
5. `pw.x` NSCF on the dense production grid
6. `kcw.x` Hamiltonian (pKI eigenvalues on the production grid)

The script is **blocking**: it polls AiiDA workchains until everything
finishes (or fails).  To run in the background:

```bash
nohup bash 00_koopmans_aiida.sh > koopmans_aiida.log 2>&1 &
tail -f koopmans_aiida.log
```

When complete, the script automatically calls `02_extract_pks.py` and prints
something like:

```
KOOPMANS_NSCF_PK = 12345   (step: .../02-wannierize/02-nscf)
KOOPMANS_HAM_PK  = 12367
```

Keep these PKs; you will need them in the next steps.

!!! note "Monitoring individual workchains"
    At any point you can inspect what is running:
    ```bash
    verdi process list -a
    verdi process status <PK>
    ```

---

## Step B — Generate the Yambo SAVE directory

Set `KOOPMANS_NSCF_PK` in `01_p2y.py` to the value printed in step A, then:

```bash
python 01_p2y.py
```

This submits a `YamboWorkflow` with `INITIALISE=True`.  The workflow runs
`pw.x` SCF + NSCF and `p2y` on the remote machine, producing the Yambo SAVE
directory.  No actual Yambo calculation is performed.

The script prints:

```
Submitted PK=<p2y_pk>  label='BSE/LDA/p2y'
Monitor with:  verdi process list -a -p <p2y_pk>
```

Wait until the workflow shows `Finished [0]` before proceeding.

---

## Step C — Extract PKs (optional)

If you need to re-run step D without re-running step A, extract the PKs
from the `step_data.pkl` directly:

```bash
python 02_extract_pks.py KOOPMANS_AIIDA/step_data.pkl
```

This prints `KOOPMANS_NSCF_PK` and `KOOPMANS_HAM_PK`.

---

## Step D — Generate ndb.QP and submit BSE

Set `P2Y_WF_PK` (from step B) and `KOOPMANS_HAM_PK` (from step A/C) in
`03_bse_submit.py`, then:

```bash
python 03_bse_submit.py
```

### What this script does

**1. Load nodes from the AiiDA database**

```python
p2y_wf       = orm.load_node(P2Y_WF_PK)
koopmans_ham = orm.load_node(KOOPMANS_HAM_PK)
```

**2. Generate ndb.QP via k2y (AiiDA calcfunction)**

```python
from k2y.aiida import generate_kcw_qp_database

kcw_db = generate_kcw_qp_database(
    yambo_retrieved,   # FolderData from the p2y YamboCalculation
    kcw_retrieved,     # FolderData from the KcwCalculation (ham)
)
```

`generate_kcw_qp_database` is an AiiDA
[`calcfunction`](https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/calculations/concepts.html#calcfunctions).
This means:

- It runs locally (no scheduler submission).
- Inputs and output (`SinglefileData` wrapping `ndb.QP`) are stored in the
  AiiDA provenance graph.
- The complete lineage from DFT → Koopmans → k2y → BSE is queryable.

**3. Build and submit the BSE YamboWorkflow**

The workflow is configured for a silicon BSE with converged parameters:

| Parameter | Value | Description |
|---|---|---|
| `BndsRnXs` | 1–100 | Screening bands in W |
| `NGsBlkXs` | 2 Ry | G-vector cutoff |
| `BSEBands`  | 4–5  | Valence + conduction bands |
| `KfnQPdb`  | `E < ./ndb.QP` | Load Koopmans QP corrections |

The `parent_folder` is set to the p2y remote folder, so the workflow reuses
the existing SAVE directory and does **not** repeat SCF + NSCF + p2y.

**4. Monitor the BSE workflow**

```bash
verdi process list -a -p <BSE_PK>
verdi process status <BSE_PK>
```

---

## Provenance graph

After all steps complete, the full provenance graph for the ndb.QP node can
be visualised with:

```bash
verdi node graph generate <kcw_db_pk>
```

The graph will show the complete lineage:

```
PwCalculation (SCF)
PwCalculation (NSCF)      KcwCalculation (ham)
      │                          │
YamboCalculation (p2y)           │
      │                          │
      └──── generate_kcw_qp_database ────┘
                    │
             SinglefileData (ndb.QP)
                    │
         YamboWorkflow (BSE)
                    │
              optical spectrum ε₂(ω)
```

---

## Expected output

The BSE YamboWorkflow stores the spectrum in its outputs.  To retrieve it:

```python
from aiida import load_profile, orm
load_profile()
bse_wf = orm.load_node(<BSE_PK>)
print(bse_wf.outputs.output_parameters.get_dict())
# {'lowest_exciton': ..., 'brightest_exciton': ..., ...}
```

The raw spectrum files (`o-bse.eps_q1_*`) are available in the retrieved
`FolderData` of the inner `YamboCalculation`.

---

## Troubleshooting

**`AssertionError` in `generate_kcw_qp_database`**
: k2y could not match all Koopmans k-points to Yambo k-points.  Ensure the
  k-point mesh in `si_ki.json` (`kpoints.grid`) matches `BSE_KMESH` in
  `01_p2y.py` and `03_bse_submit.py`.

**`p2y workflow is not finished OK`**
: Check the p2y workflow with `verdi process status <PK>` and inspect the
  logs of the failing sub-calculation.

**`Too many authentication failures` when connecting to the remote**
: See the k2y FAQ.  Use `IdentitiesOnly yes` in `~/.ssh/config` and specify
  the correct key for the computer.
