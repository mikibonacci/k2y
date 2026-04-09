# Tutorial: Silicon BSE optical spectrum with Koopmans eigenvalues

In this tutorial you will compute the optical absorption spectrum of bulk silicon
using a Bethe-Salpeter Equation (BSE) calculation in Yambo, where the quasiparticle
corrections come from a Koopmans Compliant Functional (KC) calculation rather than
from GW.

The complete input files are in `examples/silicon_scratch_to_bse/`.

## Workflow overview

```
  DFT (pw.x)
      │  SCF + NSCF on Si primitive cell
      ▼
  P2Y (p2y + yambo)
      │  Convert QE output → Yambo SAVE database
      ▼
  Koopmans (kcw.x)
      │  Compute KC eigenvalues on the k-grid
      ▼
  k2y
      │  Map KC eigenvalues → ndb.QP
      ▼
  BSE (yambo)
      │  Optical spectrum using KC quasiparticle corrections
      ▼
  optical spectrum  ε₂(ω)
```

---

## Step 0 — directory structure

The example folder is organised as follows:

```
silicon_scratch_to_bse/
├── DFT/
│   ├── scf.in       # pw.x SCF input
│   └── nscf.in      # pw.x NSCF input
├── P2Y/             # p2y conversion output (SAVE symlink)
├── KOOPMANS/        # kcw.x input/output files
├── K2Y/             # k2y output (ndb.QP)
├── YAMBO/
│   └── yambo.in     # BSE input
├── pseudo/
│   └── Si.upf
└── run.sh
```

The `run.sh` script drives all steps.  Individual steps can be selected
by passing them as arguments:

```bash
# run everything
bash run.sh

# run only K2Y and BSE (inputs already available)
bash run.sh K2Y BSE

# override binary paths without editing the script
PREFIX_YAMBO=/opt/yambo/bin bash run.sh BSE
```

---

## Step 1 — DFT: SCF and NSCF

Two `pw.x` calculations are required:

- **SCF** on a 4×4×4 k-grid to converge the charge density
- **NSCF** on a 2×2×2 grid (the one that will be used by both yambo and kcw) with 100 bands

??? note "SCF input (`DFT/scf.in`)"
    ```fortran
    &CONTROL
      calculation = 'scf'
      outdir      = './out/'
      prefix      = 'aiida'
      pseudo_dir  = '../pseudo/'
    /
    &SYSTEM
      ibrav = 0,  nat = 2,  ntyp = 1
      ecutwfc = 60,  ecutrho = 240
      force_symmorphic = .true.
    /
    &ELECTRONS
      conv_thr = 4.0d-10
    /
    ATOMIC_SPECIES
    Si  28.085  Si.upf
    ATOMIC_POSITIONS angstrom
    Si  0.000000  0.000000  0.000000
    Si -1.347153  1.347153  1.347153
    K_POINTS automatic
    4 4 4 0 0 0
    CELL_PARAMETERS angstrom
    -2.694306  0.000000  2.694306
     0.000000  2.694306  2.694306
    -2.694306  2.694306  0.000000
    ```

??? note "NSCF input (`DFT/nscf.in`)"
    ```fortran
    &CONTROL
      calculation = 'nscf'
      outdir      = './out/'
      prefix      = 'aiida'
      pseudo_dir  = './pseudo/'
    /
    &SYSTEM
      ibrav = 0,  nat = 2,  ntyp = 1
      ecutwfc = 60,  ecutrho = 240
      nbnd = 100
      force_symmorphic = .true.
    /
    &ELECTRONS
      conv_thr = 4.0d-10
      diagonalization = 'cg'
    /
    ATOMIC_SPECIES
    Si  28.085  Si.upf
    ATOMIC_POSITIONS angstrom
    Si  0.000000  0.000000  0.000000
    Si -1.347153  1.347153  1.347153
    K_POINTS automatic
    2 2 2 0 0 0
    CELL_PARAMETERS angstrom
    -2.694306  0.000000  2.694306
     0.000000  2.694306  2.694306
    -2.694306  2.694306  0.000000
    ```

Run this step with:

```bash
bash run.sh DFT
```

Output is written to `DFT/out/`.

---

## Step 2 — P2Y: create the Yambo SAVE database

`p2y` converts the Quantum ESPRESSO wavefunctions into Yambo's internal format, and
`yambo` (run without flags) initialises the remaining databases.

```bash
bash run.sh P2Y
```

After this step `P2Y/SAVE/` will contain `ns.db1`, `ndb.gops`, `ndb.kindx`, and similar files.

!!! important
    The `ns.db1` used by k2y **must** be the one from the same SAVE directory
    that will be used for the BSE calculation.  k2y reads the KS eigenvalues and
    the k-point symmetry information from that file.

---

## Step 3 — Koopmans

Run `kcw.x` to compute Koopmans Compliant eigenvalues on the same k-grid used in
the NSCF step.  The example uses the DFPT-based KC approach (`kcw-ham` mode) on
the silicon primitive cell.

```bash
bash run.sh KOOPMANS
```

The relevant output file is `KOOPMANS/Si.kcw-ham_proj.out`, which contains both
the KC (`pKI`) and KS eigenvalues on the full k-grid.

---

## Step 4 — k2y: generate the QP database

This is where k2y comes in.  It reads:

| Input | Description |
|---|---|
| `P2Y/SAVE/ns.db1` | Yambo database with KS eigenvalues and k-point symmetries |
| `KOOPMANS/Si.kcw-ham_proj.out` | KC eigenvalues from `kcw.x` |
| `KOOPMANS/Si.scf.in` | QE input used for the KC run (provides the k-point list) |

and writes `K2Y/ndb.QP`.

```bash
bash run.sh K2Y
```

which internally runs:

```bash
k2y generate \
    --ns-db1   P2Y/SAVE/ns.db1 \
    --eval     KOOPMANS/Si.kcw-ham_proj.out \
    --pwinput  KOOPMANS/Si.scf.in \
    --output   K2Y/ndb.QP
```

### What k2y does internally

1. **Load ns.db1** — reads KS eigenvalues and the irreducible k-points from Yambo.
2. **Expand to full BZ** — uses the crystal symmetry stored in ns.db1 to unfold the
   irreducible k-points to the complete Brillouin-zone grid.
3. **Load KC eigenvalues** — reads `pki_eigenvalues_on_grid` from the `kcw.x` output.
4. **Map k-points** — matches each full-BZ Yambo k-point to the corresponding KC k-point
   (with time-reversal and brute-force fallback).
5. **Compute QP correction** — for each (k, band) pair:

    $$
    \Delta E_{n\mathbf{k}} = E^{\mathrm{KI}}_{n\mathbf{k}}
    $$

    so that when Yambo applies the correction,

    $$
    E^{\mathrm{QP}}_{n\mathbf{k}} = E^{\mathrm{KS,Yambo}}_{n\mathbf{k}} + \Delta E_{n\mathbf{k}}
    = E^{\mathrm{KI}}_{n\mathbf{k}}
    $$

6. **Write ndb.QP** — creates a netCDF4 file with the correct Yambo structure.

!!! note "Version compatibility"
    k2y automatically detects the Yambo version that generated the SAVE directory
    and selects the closest bundled `ndb.QP` template.  A warning is printed if
    there is a mismatch so you know which template was used.

---

## Step 5 — BSE: optical spectrum

With `K2Y/ndb.QP` in hand, run the BSE in Yambo.  The `YAMBO/yambo.in` input already
contains the line that tells Yambo to load the QP database:

```text
KfnQPdb = 'E < ./ndb.QP'
```

The relevant BSE parameters for this silicon example are:

```text
BSEBands =  4 | 5      # valence and conduction bands included in BSE
BEnRange =  0 | 10 eV  # energy range of the spectrum
BEnSteps = 1000         # number of energy points
```

```bash
bash run.sh BSE
```

The spectrum is written to `YAMBO/bse/o-bse.eps_q1_*`.

---

## Python API

The same workflow can be driven from Python:

```python
from k2y.k2y import KcwQpDatabaseGenerator

converter = KcwQpDatabaseGenerator(ns_db1="P2Y/SAVE/ns.db1")

converter.set_koopmans_eval(path="KOOPMANS/Si.kcw-ham_proj.out")
converter.set_kpoints_from_pwinput("KOOPMANS/Si.scf.in")

converter.generate_mappings()

# optional sanity check: band gap at Γ (k_index=1), top valence = band 4
converter.verify_mappings(k_index=1, top_valence=4)

converter.generate_QP_db("K2Y/ndb.QP")
```

`verify_mappings` prints the KC and KS band gaps at the chosen k-point and raises
`AssertionError` if they are inconsistent (tolerance 1 meV).

---

## Expected output

After a successful run you will see output similar to:

```text
Reading ns.db1 from: P2Y/SAVE/ns.db1
3 kpoints expanded to 8
Loading Koopmans eigenvalues from: KOOPMANS/Si.kcw-ham_proj.out
Loading k-points from: KOOPMANS/Si.scf.in
KcwQpDatabaseGenerator Summary
========================================
ns.db1 loaded: True
  K-points in ns.db1: 1
Koopmans eigenvalues: (1, 160)
KCW k-points: 8 (crystal)
Mappings generated: No
Generating mappings...
Writing QP database to: K2Y/ndb.QP
...done.
Done.
```

The generated `ndb.QP` will have:

| Variable | Shape | Description |
|---|---|---|
| `QP_kpts` | (3, 8) | Full-BZ k-points |
| `QP_table` | (3, 160) | Band/k-point index table |
| `QP_E` | (160, 2) | KC quasiparticle energies (Re, Im) |
| `QP_Eo` | (160,) | Reference KS energies |
