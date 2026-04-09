# Quick Start

## Prerequisites

k2y requires:

- Python ≥ 3.8
- A working installation of [Quantum ESPRESSO](https://www.quantum-espresso.org/) (for `pw.x` and `kcw.x`)
- A working installation of [Yambo](https://www.yambo-code.eu/) (for `p2y` and `yambo`)

The following Python packages are **not on PyPI** and must be installed from source before k2y:

| Package | Source | Notes |
|---|---|---|
| `ase-koopmans` | `https://github.com/mikibonacci/ase_koopmans` branch `add/pKI_diag` | Koopmans fork of ASE |
| `yambopy` | `https://github.com/yambo-code/yambopy` branch `aiida-parsing` | Python interface to Yambo |

```bash
pip install git+https://github.com/mikibonacci/ase_koopmans.git@add/pKI_diag#egg=ase-koopmans
pip install git+https://github.com/yambo-code/yambopy.git@aiida-parsing#egg=yambopy
```

## Installation

### From source (recommended)

```bash
git clone https://github.com/mikibonacci/k2y.git
cd k2y
pip install -e .
```

### With optional AiiDA support

```bash
pip install -e ".[aiida]"
```

### Verify the installation

```bash
k2y --help
```

You should see:

```text
Usage: k2y [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  generate   Generate a Yambo ndb.QP from Koopmans eigenvalues.
  kpoints    Print the k-points card needed for a Koopmans interpolation run.
```

## Minimal example

If you already have a Yambo `SAVE` directory and a `kcw.x` output file, generating
the QP database is a single command:

```bash
k2y generate \
    --ns-db1   SAVE/ns.db1 \
    --eval     Si.kcw-ham_proj.out \
    --pwinput  Si.scf.in \
    --output   ndb.QP
```

Then place `ndb.QP` in your Yambo run directory and point to it in `yambo.in`:

```text
KfnQPdb = 'E < ./ndb.QP'
```

See the [Silicon BSE tutorial](../tutorials/silicon_bse.md) for a complete worked example.
