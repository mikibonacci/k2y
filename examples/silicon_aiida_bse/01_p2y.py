#!/usr/bin/env python
"""
01_p2y.py – submit a YamboWorkflow to generate the Yambo SAVE directory
=========================================================================

Step B of the AiiDA-based BSE@Koopmans workflow.

This script submits a YamboWorkflow with INITIALISE=True, which runs:
  QE scf → QE nscf → p2y → yambo (init)

No Yambo BSE calculation is performed.  The sole purpose is to produce the
Yambo SAVE directory on the remote machine so that the BSE step (03_bse_submit.py)
can reuse it without repeating the expensive scf/nscf.

After running:
  1. Wait for the submitted workflow to finish  (verdi process list -a -p <PK>)
  2. Set P2Y_WF_PK in 03_bse_submit.py to the PK printed at the end.

Prerequisites
-------------
  - Step A (00_koopmans_aiida.sh) must have finished.
  - Set KOOPMANS_NSCF_PK below (from 02_extract_pks.py output).
  - AiiDA codes and pseudopotential family must be registered.

Usage
-----
  python 01_p2y.py
"""

from aiida import load_profile, orm
from aiida.engine import submit
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType

load_profile()

# ---------------------------------------------------------------------------
# User settings – adapt for your system / cluster
# ---------------------------------------------------------------------------

# PK of the PwBaseWorkChain from the Koopmans NSCF step.
# Obtain it by running 02_extract_pks.py against the step_data.pkl produced
# by 00_koopmans_aiida.sh.
KOOPMANS_NSCF_PK = None   # <-- set this (printed by 02_extract_pks.py)

# k-point mesh for the p2y / BSE run (must match the Koopmans NSCF mesh)
BSE_KMESH = [2, 2, 2]

# AiiDA code labels (as registered via `verdi code`)
PW_CODE            = "pw@localhost"          # <-- adapt
PREPROCESSING_CODE = "p2y@localhost"         # <-- adapt
YAMBO_CODE         = "yambo@localhost"       # <-- adapt

# Pseudo-potential family
PSEUDO_FAMILY = "PseudoDojo/0.4/LDA/SR/standard/upf"

# Scheduler resources
METADATA_YAMBO = {"options": {
    "max_wallclock_seconds": 86400,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 12,
        "num_cores_per_mpiproc": 1,
    },
    "custom_scheduler_commands": "export OMP_NUM_THREADS=1",
}}

METADATA_SCF = {"options": {
    "max_wallclock_seconds": 86400,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 12,
        "num_cores_per_mpiproc": 1,
    },
    "custom_scheduler_commands": "export OMP_NUM_THREADS=1",
}}

METADATA_NSCF = {"options": {
    "max_wallclock_seconds": 86400,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 12,
        "num_cores_per_mpiproc": 1,
    },
    "custom_scheduler_commands": "export OMP_NUM_THREADS=1",
}}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
if KOOPMANS_NSCF_PK is None:
    raise ValueError(
        "KOOPMANS_NSCF_PK is not set.  Run 02_extract_pks.py against the "
        "step_data.pkl produced by 00_koopmans_aiida.sh and set the value here."
    )

# ---------------------------------------------------------------------------
# Load structure and ecutwfc from Koopmans NSCF node
# ---------------------------------------------------------------------------
print(f"Loading Koopmans NSCF node PK={KOOPMANS_NSCF_PK}")
koopmans_nscf = orm.load_node(KOOPMANS_NSCF_PK)
structure = koopmans_nscf.inputs.pw.structure
ecutwfc   = koopmans_nscf.inputs.pw.parameters["SYSTEM"]["ecutwfc"]
print(f"  structure : {structure.get_formula()}")
print(f"  ecutwfc   : {ecutwfc} Ry")

kpoints = orm.KpointsData()
kpoints.set_kpoints_mesh(BSE_KMESH)

# ---------------------------------------------------------------------------
# Build the YamboWorkflow for p2y only (INITIALISE=True)
# ---------------------------------------------------------------------------
YamboWorkflow = WorkflowFactory("yambo.yambo.yambowf")

overrides = {
    "yres": {
        "yambo": {
            "metadata": METADATA_YAMBO,
        },
    },
    "nscf": {
        "pw": {
            "metadata": METADATA_NSCF,
            "parameters": {
                "SYSTEM":    {"ecutwfc": ecutwfc},
                "ELECTRONS": {"diagonalization": "cg"},
            },
        },
    },
    "scf": {
        "pw": {
            "metadata": METADATA_SCF,
            "parameters": {
                "SYSTEM": {"ecutwfc": ecutwfc},
            },
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

# INITIALISE=True: stop after scf + nscf + p2y; don't run any Yambo step
builder.yres.yambo.settings = orm.Dict(dict={"INITIALISE": True})

# Override k-points
builder.nscf.kpoints = kpoints

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
print("\nSubmitting YamboWorkflow for p2y (INITIALISE=True) ...")
node = submit(builder)
node.label = "BSE/LDA/p2y"
print(f"  Submitted PK={node.pk}  label='{node.label}'")
print(f"  Monitor with:  verdi process list -a -p {node.pk}")
print()
print("Once finished OK, set P2Y_WF_PK in 03_bse_submit.py:")
print(f"  P2Y_WF_PK = {node.pk}")
