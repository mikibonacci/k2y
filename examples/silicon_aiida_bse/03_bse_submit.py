#!/usr/bin/env python
"""
03_bse_submit.py – generate ndb.QP via k2y and submit BSE YamboWorkflow
=========================================================================

Step D of the AiiDA-based BSE@Koopmans workflow.  Prerequisites:

  Step A  00_koopmans_aiida.sh  finished OK  → step_data.pkl produced
  Step B  01_p2y.py             finished OK  → P2Y_WF_PK known
  Step C  02_extract_pks.py     run          → KOOPMANS_HAM_PK known

This script:
  1. Loads the completed p2y YamboWorkflow and extracts its outputs.
  2. Loads the KcwCalculation (ham) node from the Koopmans pKI run.
  3. Generates the Yambo ndb.QP database via k2y.aiida.generate_kcw_qp_database
     (an AiiDA calcfunction – provenance is tracked automatically).
  4. Builds and submits a YamboWorkflow for BSE, reusing the p2y remote folder
     so that scf + nscf + p2y are not repeated.

Usage
-----
  python 03_bse_submit.py

Set P2Y_WF_PK and KOOPMANS_HAM_PK at the top of this file.
"""

from aiida import load_profile, orm
from aiida.engine import submit
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType

from k2y.aiida import generate_kcw_qp_database

load_profile()

# ---------------------------------------------------------------------------
# User settings – adapt for your system / cluster
# ---------------------------------------------------------------------------

# PK of the completed p2y YamboWorkflow (output of 01_p2y.py)
P2Y_WF_PK = None          # <-- set this after 01_p2y.py finishes OK

# PK of the KcwCalculation (ham) from the Koopmans pKI run
# (from 02_extract_pks.py)
KOOPMANS_HAM_PK = None    # <-- set this

# Optional: pre-existing ndb.QP node to use as template for k2y (copies
# dimensions/attributes).  Set to None to use the bundled k2y template.
QP_TEMPLATE_PK = None

# AiiDA code labels (as registered via `verdi code`)
PW_CODE            = "pw@localhost"          # <-- adapt
PREPROCESSING_CODE = "p2y@localhost"         # <-- adapt
YAMBO_CODE         = "yambo@localhost"       # <-- adapt

# Pseudo-potential family
PSEUDO_FAMILY = "PseudoDojo/0.4/LDA/SR/standard/upf"

# k-point mesh for the BSE calculation (must match the NSCF mesh)
BSE_KMESH = [2, 2, 2]

# BSE parameters (converged values for silicon)
NB_BANDS  = 100    # screening bands in W (BndsRnXs)
GCUT_RY   = 2      # G-vector cutoff in Ry (NGsBlkXs, BSENGBlk)
BSE_BANDS = [4, 5] # [first valence, last conduction] band window for BSE

# Scheduler resources
METADATA_YAMBO = {"options": {
    "max_wallclock_seconds": 86400,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 4,
        "num_cores_per_mpiproc": 1,
    },
    "custom_scheduler_commands": "export OMP_NUM_THREADS=1",
}}

METADATA_SCF = {"options": {
    "max_wallclock_seconds": 86400,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 4,
        "num_cores_per_mpiproc": 1,
    },
    "custom_scheduler_commands": "export OMP_NUM_THREADS=1",
}}

METADATA_NSCF = {"options": {
    "max_wallclock_seconds": 86400,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 4,
        "num_cores_per_mpiproc": 1,
    },
    "custom_scheduler_commands": "export OMP_NUM_THREADS=1",
}}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
if P2Y_WF_PK is None:
    raise ValueError(
        "P2Y_WF_PK is not set.  Run 01_p2y.py first, wait for it to finish "
        "successfully, then set P2Y_WF_PK to the PK printed by that script."
    )
if KOOPMANS_HAM_PK is None:
    raise ValueError(
        "KOOPMANS_HAM_PK is not set.  Run 02_extract_pks.py and copy the "
        "KOOPMANS_HAM_PK value here."
    )

# ---------------------------------------------------------------------------
# Step 1 – load nodes
# ---------------------------------------------------------------------------
print(f"Loading p2y YamboWorkflow node PK={P2Y_WF_PK}")
p2y_wf = orm.load_node(P2Y_WF_PK)
assert p2y_wf.is_finished_ok, (
    f"p2y workflow PK={P2Y_WF_PK} is not finished OK "
    f"(exit_status={p2y_wf.exit_status}).  Cannot continue."
)

# FolderData produced by the inner YamboCalculation – required by k2y
yambo_retrieved = p2y_wf.outputs.retrieved
print(f"  yambo_retrieved PK: {yambo_retrieved.pk}")

# RemoteData for the p2y folder – reused to skip scf + nscf + p2y in the BSE step
p2y_remote = yambo_retrieved.creator.outputs.remote_folder
print(f"  p2y remote folder PK: {p2y_remote.pk}")

# Structure and ecutwfc – taken from the p2y workflow inputs for consistency
structure = p2y_wf.inputs.nscf__pw__structure
ecutwfc   = p2y_wf.inputs.nscf__pw__parameters["SYSTEM"].get("ecutwfc", None)
print(f"  structure : {structure.get_formula()}")

print(f"\nLoading Koopmans HAM node PK={KOOPMANS_HAM_PK}")
koopmans_ham  = orm.load_node(KOOPMANS_HAM_PK)
kcw_retrieved = koopmans_ham.outputs.retrieved
print(f"  kcw_retrieved PK: {kcw_retrieved.pk}")

# ---------------------------------------------------------------------------
# Step 2 – generate Yambo ndb.QP via k2y (AiiDA calcfunction)
# ---------------------------------------------------------------------------
print("\nGenerating Koopmans ndb.QP database via k2y ...")
qp_template = orm.load_node(QP_TEMPLATE_PK) if QP_TEMPLATE_PK else None

kcw_db = generate_kcw_qp_database(
    yambo_retrieved,
    kcw_retrieved,
    QP_template_node=qp_template,
)
print(f"  ndb.QP SinglefileData PK: {kcw_db.pk}")

# ---------------------------------------------------------------------------
# Step 3 – build the YamboWorkflow for BSE
# ---------------------------------------------------------------------------
print("\nBuilding YamboWorkflow builder for BSE ...")

YamboWorkflow = WorkflowFactory("yambo.yambo.yambowf")

kpoints = orm.KpointsData()
kpoints.set_kpoints_mesh(BSE_KMESH)

overrides_pw = {}
if ecutwfc:
    overrides_pw = {"SYSTEM": {"ecutwfc": ecutwfc}}

overrides_yambo = {
    "nscf": {
        "pw": {
            "parameters": {
                **overrides_pw,
                "ELECTRONS": {"diagonalization": "cg"},
            },
        },
    },
    "scf": {
        "pw": {
            "parameters": overrides_pw,
        },
    },
    "yambo": {
        "parameters": {
            "arguments": [
                "rim_cut",
                "WRbsWF",
                "NLCC",   # required for PseudoDojo pseudopotentials
            ],
            "variables": {
                "BndsRnXs": [[1, NB_BANDS], ""],
                "NGsBlkXs": [GCUT_RY, "Ry"],
                "BSENGBlk": [GCUT_RY, "Ry"],
                "BSEBands": [BSE_BANDS, ""],
                "LongDrXs": [[1.0, 1.0, 1.0], ""],
                "BLongDir": [[1.0, 1.0, 1.0], ""],
                "BEnRange": [[0, 10], "eV"],
                "BEnSteps": [1000, ""],
                "BDmRange": [[0.1, 0.1], "eV"],
                # Load Koopmans QP corrections
                "KfnQPdb": "E < ./ndb.QP",
                "BS_ROLEs":"k eh t",
                "BS_CPU":"2 2 1",
            },
        },
        "metadata": METADATA_YAMBO,
    },
}

overrides = {
    "yres": overrides_yambo,
    "nscf": {
        "pw": {
            "metadata":   METADATA_NSCF,
            "parameters": overrides_yambo["nscf"]["pw"]["parameters"],
        },
    },
    "scf": {
        "pw": {
            "metadata":   METADATA_SCF,
            "parameters": overrides_yambo["scf"]["pw"]["parameters"],
        },
    },
}

builder = YamboWorkflow.get_builder_from_protocol(
    pw_code            = PW_CODE,
    preprocessing_code = PREPROCESSING_CODE,
    code               = YAMBO_CODE,
    protocol           = "moderate",
    protocol_qe        = "moderate",
    structure          = structure,
    pseudo_family      = PSEUDO_FAMILY,
    overrides          = overrides,
    electronic_type    = ElectronicType.INSULATOR,
    calc_type          = "bse",
)

# Request parsing of lowest and brightest exciton energies
builder.additional_parsing = orm.List(list=["lowest_exciton", "brightest_exciton"])

# Override k-points with the BSE mesh
builder.nscf.kpoints = kpoints

# Attach the Koopmans QP database generated by k2y
builder.yres.yambo.QP_corrections = kcw_db

# Remove GW-specific keys that get_builder_from_protocol may inject
bse_params = builder.yres.yambo.parameters.get_dict()
for key in ("GTermKinds", "GbndRnge", "FFTGvecs"):
    bse_params["variables"].pop(key, None)
builder.yres.yambo.parameters = orm.Dict(dict=bse_params)

# Reuse the p2y remote folder – skips scf + nscf + p2y inside the workflow
if not p2y_remote.is_empty:
    builder.parent_folder = p2y_remote
    print(f"  Reusing p2y remote folder PK={p2y_remote.pk}")
else:
    builder.yres.yambo.settings = orm.Dict(dict={"INITIALISE": True})
    print("  WARNING: p2y remote folder is empty – workflow will re-run p2y.")

# ---------------------------------------------------------------------------
# Step 4 – submit
# ---------------------------------------------------------------------------
print("\nSubmitting YamboWorkflow for BSE@pKI ...")
node = submit(builder)
node.label = "BSE/LDA/pKI"
print(f"  Submitted PK={node.pk}  label='{node.label}'")
print(f"  Monitor with:  verdi process list -a -p {node.pk}")
