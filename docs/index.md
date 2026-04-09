# k2y — Koopmans to Yambo

**k2y** converts [Koopmans Functionals](https://koopmans-functionals.org/) eigenvalues from
`kcw.x` (Quantum ESPRESSO) into Yambo-compatible quasiparticle databases (`ndb.QP`), enabling
BSE optical-spectrum calculations without the cost of a GW run.

[![PyPI version](https://badge.fury.io/py/k2y.svg)](https://badge.fury.io/py/k2y)
[![Build Status](https://github.com/mikibonacci/k2y/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/mikibonacci/k2y/actions)
[![Docs status](https://github.com/mikibonacci/k2y/actions/workflows/docs.yml/badge.svg)](https://mikibonacci.github.io/k2y)

---

<div class="grid cards" markdown>

-   **Getting Started**

    ---

    Install k2y and run your first conversion in minutes.

    [To the quick start guide](installation_quickstart/quickstart.md)

-   **Tutorials**

    ---

    Step-by-step walkthroughs, starting from a Silicon BSE example.

    [To the tutorials](tutorials/index.md)

</div>

---

## What does k2y do?

A standard many-body perturbation theory workflow for optical properties looks like:

```
DFT (pw.x) → GW (yambo) → BSE (yambo)
```

Koopmans Compliant (KC) functionals provide quasiparticle-quality band structures at DFT cost.
**k2y** bridges the two codes so you can replace the expensive GW step:

```
DFT (pw.x) → Koopmans (kcw.x) → k2y → BSE (yambo)
```

Internally, k2y:

1. Reads the Yambo `ns.db1` database (KS eigenvalues, k-point grid, symmetries)
2. Loads KC eigenvalues from the `kcw.x` output
3. Maps k-points between the two grids, expanding from the irreducible BZ to the full one
4. Writes a `ndb.QP` file that Yambo reads as QP corrections during a BSE run

## How to cite

If you use k2y in your research, please cite:

> M. Bonacci *et al.*, *k2y: bridging Koopmans Functionals and Yambo for optical spectra*, (in preparation).

## Acknowledgements

Development of k2y was supported by the [MaX — Materials Design at the Exascale](http://www.max-centre.eu/) Centre of Excellence.
