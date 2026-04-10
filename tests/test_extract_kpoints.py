"""
Tests for k2y.extract_kpoints.extract_kpoints_from_pwin.
"""
import textwrap

import numpy as np
import pytest

from k2y.extract_kpoints import extract_kpoints_from_pwin


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_pwin(tmp_path, kpoints_block: str) -> str:
    """Write a minimal pw.x input containing kpoints_block and return the path."""
    content = textwrap.dedent("""\
        &CONTROL
          calculation = 'nscf'
        /
        &SYSTEM
          ibrav = 0, nat = 2, ntyp = 1
          ecutwfc = 20.0
        /
        &ELECTRONS
        /
        ATOMIC_SPECIES
        Si 28.085 Si.upf
        ATOMIC_POSITIONS crystal
        Si 0.00 0.00 0.00
        Si 0.25 0.25 0.25
    """) + "\n" + textwrap.dedent(kpoints_block)
    p = tmp_path / "test.in"
    p.write_text(content)
    return str(p)


# ---------------------------------------------------------------------------
# Tests against the real Si.scf.in fixture
# ---------------------------------------------------------------------------

class TestRealSiScfIn:

    def test_returns_ndarray(self, scf_in_path):
        kpts, _ = extract_kpoints_from_pwin(str(scf_in_path))
        assert isinstance(kpts, np.ndarray)

    def test_kpoints_count(self, scf_in_path):
        """Si.scf.in contains 8 explicit k-points."""
        kpts, _ = extract_kpoints_from_pwin(str(scf_in_path))
        assert kpts.shape[0] == 8

    def test_kpoints_shape(self, scf_in_path):
        kpts, _ = extract_kpoints_from_pwin(str(scf_in_path))
        assert kpts.ndim == 2
        assert kpts.shape[1] == 3

    def test_kpoints_type_crystal(self, scf_in_path):
        _, ktype = extract_kpoints_from_pwin(str(scf_in_path))
        assert ktype.lower() == "crystal"

    def test_gamma_point_present(self, scf_in_path):
        """The 2×2×2 grid for Si must include Γ = (0, 0, 0)."""
        kpts, _ = extract_kpoints_from_pwin(str(scf_in_path))
        gamma = np.array([0.0, 0.0, 0.0])
        dists = np.linalg.norm(kpts - gamma, axis=1)
        assert dists.min() < 1e-6

    def test_kpoints_in_unit_cell(self, scf_in_path):
        """All k-point coordinates must be in [0, 1) for a crystal-type grid."""
        kpts, _ = extract_kpoints_from_pwin(str(scf_in_path))
        assert np.all(kpts >= -1e-6)
        assert np.all(kpts <= 1.0 + 1e-6)


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        extract_kpoints_from_pwin("/nonexistent/input.in")


def test_no_kpoints_card(tmp_path):
    content = textwrap.dedent("""\
        &CONTROL
        /
        &SYSTEM
          ibrav = 0, nat = 1, ntyp = 1
        /
    """)
    p = tmp_path / "no_kpts.in"
    p.write_text(content)
    with pytest.raises(ValueError, match="K_POINTS"):
        extract_kpoints_from_pwin(str(p))


# ---------------------------------------------------------------------------
# Synthetic inputs
# ---------------------------------------------------------------------------

def test_tpiba_type(tmp_path):
    block = """\
        K_POINTS tpiba
        2
        0.0  0.0  0.0  0.5
        0.5  0.0  0.0  0.5
    """
    p = _write_pwin(tmp_path, block)
    kpts, ktype = extract_kpoints_from_pwin(p)
    assert ktype.lower() == "tpiba"
    assert kpts.shape == (2, 3)
    np.testing.assert_allclose(kpts[0], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(kpts[1], [0.5, 0.0, 0.0])


def test_weights_not_included_in_output(tmp_path):
    """Third column (weights) should be excluded; only first 3 coordinates kept."""
    block = """\
        K_POINTS crystal
        3
        0.0  0.0  0.0  0.125
        0.5  0.0  0.0  0.125
        0.5  0.5  0.0  0.125
    """
    p = _write_pwin(tmp_path, block)
    kpts, _ = extract_kpoints_from_pwin(p)
    assert kpts.shape == (3, 3)
    np.testing.assert_allclose(kpts[2], [0.5, 0.5, 0.0])
