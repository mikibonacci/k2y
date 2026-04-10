"""
Tests for the k2y command-line interface (CLI).

Uses Click's test runner so that no subprocess is needed and coverage is
properly tracked.
"""
import netCDF4 as nc
import numpy as np
import pytest
from click.testing import CliRunner

from k2y.cli import main

N_FULL_KPTS = 64
N_BANDS     = 20
N_STATES    = N_FULL_KPTS * N_BANDS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@pytest.fixture()
def runner():
    return CliRunner()


# ---------------------------------------------------------------------------
# Top-level help / version
# ---------------------------------------------------------------------------

class TestTopLevel:

    def test_help_exits_zero(self, runner):
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0

    def test_help_mentions_generate(self, runner):
        result = runner.invoke(main, ["--help"])
        assert "generate" in result.output

    def test_version_exits_zero(self, runner):
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0

    def test_generate_help_exits_zero(self, runner):
        result = runner.invoke(main, ["generate", "--help"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# generate – file mode
# ---------------------------------------------------------------------------

class TestGenerateFileMode:

    def test_creates_output_file(
        self, runner, ns_db1_path, kcw_out_path, scf_in_path, tmp_path
    ):
        out = str(tmp_path / "ndb.QP")
        result = runner.invoke(main, [
            "generate",
            "--ns-db1",   str(ns_db1_path),
            "--eval",     str(kcw_out_path),
            "--pwinput",  str(scf_in_path),
            "--output",   out,
        ])
        assert result.exit_code == 0, result.output
        import os; assert os.path.exists(out)

    def test_output_has_correct_shapes(
        self, runner, ns_db1_path, kcw_out_path, scf_in_path, tmp_path
    ):
        out = str(tmp_path / "ndb.QP")
        runner.invoke(main, [
            "generate",
            "--ns-db1",   str(ns_db1_path),
            "--eval",     str(kcw_out_path),
            "--pwinput",  str(scf_in_path),
            "--output",   out,
        ])
        with nc.Dataset(out) as ds:
            assert ds.variables["QP_E"][:].shape    == (N_STATES, 2)
            assert ds.variables["QP_Eo"][:].shape   == (N_STATES,)
            assert ds.variables["QP_kpts"][:].shape == (3, N_FULL_KPTS)

    def test_output_pars_correct(
        self, runner, ns_db1_path, kcw_out_path, scf_in_path, tmp_path
    ):
        out = str(tmp_path / "ndb.QP")
        runner.invoke(main, [
            "generate",
            "--ns-db1",   str(ns_db1_path),
            "--eval",     str(kcw_out_path),
            "--pwinput",  str(scf_in_path),
            "--output",   out,
        ])
        with nc.Dataset(out) as ds:
            pars = np.array(ds.variables["PARS"][:]).flatten()
        assert int(pars[0]) == N_BANDS
        assert int(pars[1]) == N_FULL_KPTS

    def test_with_verify_flags(
        self, runner, ns_db1_path, kcw_out_path, scf_in_path, tmp_path
    ):
        """verify_mappings should not cause the CLI to crash."""
        out = str(tmp_path / "ndb.QP")
        result = runner.invoke(main, [
            "generate",
            "--ns-db1",   str(ns_db1_path),
            "--eval",     str(kcw_out_path),
            "--pwinput",  str(scf_in_path),
            "--output",   out,
            "--verify-k",  "1",
            "--verify-tv", "4",
        ])
        assert result.exit_code == 0, result.output

    def test_default_output_name(
        self, runner, ns_db1_path, kcw_out_path, scf_in_path, tmp_path
    ):
        """When --output is omitted, ndb.QP is written in the current directory."""
        import os
        with runner.isolated_filesystem(temp_dir=str(tmp_path)):
            result = runner.invoke(main, [
                "generate",
                "--ns-db1",  str(ns_db1_path),
                "--eval",    str(kcw_out_path),
                "--pwinput", str(scf_in_path),
            ])
            assert result.exit_code == 0, result.output
            assert os.path.exists("ndb.QP")


# ---------------------------------------------------------------------------
# generate – missing argument errors
# ---------------------------------------------------------------------------

class TestGenerateMissingArgs:

    def test_missing_ns_db1_exits_nonzero(
        self, runner, kcw_out_path, scf_in_path, tmp_path
    ):
        result = runner.invoke(main, [
            "generate",
            "--eval",    str(kcw_out_path),
            "--pwinput", str(scf_in_path),
            "--output",  str(tmp_path / "ndb.QP"),
        ])
        assert result.exit_code != 0

    def test_missing_eval_exits_nonzero(
        self, runner, ns_db1_path, scf_in_path, tmp_path
    ):
        result = runner.invoke(main, [
            "generate",
            "--ns-db1",  str(ns_db1_path),
            "--pwinput", str(scf_in_path),
            "--output",  str(tmp_path / "ndb.QP"),
        ])
        assert result.exit_code != 0

    def test_missing_pwinput_exits_nonzero(
        self, runner, ns_db1_path, kcw_out_path, tmp_path
    ):
        result = runner.invoke(main, [
            "generate",
            "--ns-db1", str(ns_db1_path),
            "--eval",   str(kcw_out_path),
            "--output", str(tmp_path / "ndb.QP"),
        ])
        assert result.exit_code != 0

    def test_nonexistent_ns_db1_exits_nonzero(
        self, runner, kcw_out_path, scf_in_path, tmp_path
    ):
        result = runner.invoke(main, [
            "generate",
            "--ns-db1",  "/nonexistent/ns.db1",
            "--eval",    str(kcw_out_path),
            "--pwinput", str(scf_in_path),
            "--output",  str(tmp_path / "ndb.QP"),
        ])
        assert result.exit_code != 0
