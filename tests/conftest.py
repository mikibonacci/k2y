"""
Shared pytest fixtures for the k2y test suite.
"""
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Return the absolute path to the tests/data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def ns_db1_path(data_dir) -> Path:
    return data_dir / "ns.db1"


@pytest.fixture(scope="session")
def kcw_out_path(data_dir) -> Path:
    return data_dir / "Si.kcw-ham_proj.out"


@pytest.fixture(scope="session")
def scf_in_path(data_dir) -> Path:
    return data_dir / "Si.scf.in"


@pytest.fixture(scope="session")
def ref_ndb_qp_path(data_dir) -> Path:
    return data_dir / "ref_ndb.QP"


# ---------------------------------------------------------------------------
# Pre-built generator (with mappings) – session-scope to avoid re-reading
# large files in every test
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def generator(ns_db1_path, kcw_out_path, scf_in_path):
    """Fully initialised KcwQpDatabaseGenerator with mappings already generated."""
    import warnings
    from k2y.k2y import KcwQpDatabaseGenerator

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_koopmans_eval(path=str(kcw_out_path))
        gen.set_kpoints_from_pwinput(str(scf_in_path))
        gen.generate_mappings()
    return gen


@pytest.fixture()
def tmp_qp_path(tmp_path) -> Path:
    """A temporary file path for the output ndb.QP."""
    return tmp_path / "ndb.QP"
