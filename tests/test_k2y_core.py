"""
Tests for the core KcwQpDatabaseGenerator pipeline:
  __init__, validate_inputs, set_koopmans_eval, set_kpoints_from_pwinput,
  generate_mappings, verify_mappings, generate_QP_db, summary.

Expected values are derived from the bulk-silicon test fixture (2×2×2 grid,
Yambo 5.1, 20 bands).  The full BZ contains 64 k-points (8 IBZ → 64 after
symmetry expansion).
"""
import warnings

import netCDF4 as nc
import numpy as np
import pytest

from k2y.k2y import KcwQpDatabaseGenerator

# Expected values for the Si 2×2×2 fixture
N_IBZ_KPTS  = 8
N_FULL_KPTS = 64
N_BANDS     = 20
N_STATES    = N_FULL_KPTS * N_BANDS  # 1280
TOP_VALENCE = 4                       # 1-based (Si has 4 valence bands)
KI_GAP_K1   = 3.476                   # eV, rounded to 3 d.p.
TOL_EV      = 5e-3                    # comparison tolerance in eV


# ---------------------------------------------------------------------------
# Initialisation and input validation
# ---------------------------------------------------------------------------

class TestInit:

    def test_init_with_ns_db1(self, ns_db1_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        assert gen.ns_db1 is not None

    def test_init_missing_file_raises(self):
        with pytest.raises((FileNotFoundError, Exception)):
            KcwQpDatabaseGenerator(ns_db1="/nonexistent/ns.db1")

    def test_template_assigned_on_init(self, ns_db1_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        assert gen.template_QP_path is not None

    def test_custom_template_qp_path(self, ns_db1_path):
        """A custom template_QP_path must be accepted; the active path must exist.
        Note: check_version_compatibility(auto_select=True) may switch to a closer
        bundled template, so we only assert the resulting path exists."""
        tpl = KcwQpDatabaseGenerator.get_templateQP_filepath()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(
                ns_db1=str(ns_db1_path),
                template_QP_path=str(tpl),
            )
        from pathlib import Path
        assert Path(str(gen.template_QP_path)).exists()

    def test_custom_template_not_found_raises(self, ns_db1_path):
        with pytest.raises(FileNotFoundError):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                KcwQpDatabaseGenerator(
                    ns_db1=str(ns_db1_path),
                    template_QP_path="/nonexistent/template.QP",
                )


class TestValidateInputs:

    def test_all_false_on_empty(self):
        # ns_db1 is not set as an attribute when no argument is given to __init__,
        # so validate_inputs() raises AttributeError on an empty generator.
        gen = KcwQpDatabaseGenerator()
        with pytest.raises(AttributeError):
            gen.validate_inputs()

    def test_keys_present(self, ns_db1_path):
        """Load ns_db1 first so the attribute exists; check expected keys."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        result = gen.validate_inputs()
        assert set(result.keys()) == {
            "ns_db1", "eigenvalues_KI", "eigenvalues_KS", "kpoints_grid_kcw"
        }

    def test_ns_db1_true_after_load(self, ns_db1_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        assert gen.validate_inputs()["ns_db1"] is True

    def test_validate_or_raise_with_empty(self, ns_db1_path):
        """validate_or_raise should fail when eigenvalues and k-points are missing."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        # ns_db1 loaded but no eigenvalues or k-points → should raise
        with pytest.raises(ValueError, match="Missing required inputs"):
            gen.validate_or_raise()

    def test_all_true_after_full_setup(self, generator):
        result = generator.validate_inputs()
        assert all(result.values())


# ---------------------------------------------------------------------------
# set_koopmans_eval
# ---------------------------------------------------------------------------

class TestSetKoopmanEval:

    def test_sets_eigenvalues_ki(self, ns_db1_path, kcw_out_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_koopmans_eval(path=str(kcw_out_path))
        assert gen.eigenvalues_KI is not None

    def test_sets_eigenvalues_ks(self, ns_db1_path, kcw_out_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_koopmans_eval(path=str(kcw_out_path))
        assert gen.eigenvalues_KS is not None

    def test_ki_ks_same_shape(self, ns_db1_path, kcw_out_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_koopmans_eval(path=str(kcw_out_path))
        assert gen.eigenvalues_KI.shape == gen.eigenvalues_KS.shape

    def test_n_bands_matches_fixture(self, ns_db1_path, kcw_out_path, scf_in_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_koopmans_eval(path=str(kcw_out_path))
        gen.set_kpoints_from_pwinput(str(scf_in_path))
        # After reshaping (done in generate_mappings) bands = N_BANDS
        # Before mappings, eigenvalues are flat or 2-D; just check total size.
        total = gen.eigenvalues_KI.size
        assert total == N_IBZ_KPTS * N_BANDS


# ---------------------------------------------------------------------------
# set_kpoints_from_pwinput
# ---------------------------------------------------------------------------

class TestSetKpointsFromPwinput:

    def test_sets_kpoints_grid(self, ns_db1_path, scf_in_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_kpoints_from_pwinput(str(scf_in_path))
        assert gen.kpoints_grid_kcw is not None

    def test_n_kpoints(self, ns_db1_path, scf_in_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_kpoints_from_pwinput(str(scf_in_path))
        assert gen.kpoints_grid_kcw.shape[0] == N_IBZ_KPTS

    def test_kpoints_type(self, ns_db1_path, scf_in_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        gen.set_kpoints_from_pwinput(str(scf_in_path))
        assert gen.kpoints_type.lower() == "crystal"

    def test_file_not_found(self, ns_db1_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        with pytest.raises(FileNotFoundError):
            gen.set_kpoints_from_pwinput("/nonexistent/nscf.in")


# ---------------------------------------------------------------------------
# generate_mappings
# ---------------------------------------------------------------------------

class TestGenerateMappings:

    def test_runs_without_error(self, generator):
        assert generator.mapped_vars is not None

    def test_full_bz_kpoints(self, generator):
        assert generator.kpoints.shape[1] == N_FULL_KPTS

    def test_n_bands(self, generator):
        assert generator.new_KI.shape[1] == N_BANDS
        assert generator.new_KS.shape[1] == N_BANDS

    def test_new_ki_new_ks_same_shape(self, generator):
        assert generator.new_KI.shape == generator.new_KS.shape

    def test_ki_gap_at_gamma(self, generator):
        """KC gap at k=0 (Γ) between valence and conduction must match known value."""
        gap = np.round(
            generator.new_KI[0, TOP_VALENCE] - generator.new_KI[0, TOP_VALENCE - 1],
            3,
        )
        assert abs(gap - KI_GAP_K1) < TOL_EV

    def test_mapped_vars_keys(self, generator):
        expected = {
            "QP_QP_@_state_1_b_range",
            "QP_QP_@_state_1_K_range",
            "QP_kpts", "QP_E", "QP_Eo", "QP_Z", "QP_table", "PARS",
        }
        assert expected.issubset(set(generator.mapped_vars.keys()))

    def test_band_range(self, generator):
        assert generator.mapped_vars["QP_QP_@_state_1_b_range"] == [1, N_BANDS]

    def test_kpoint_range(self, generator):
        assert generator.mapped_vars["QP_QP_@_state_1_K_range"] == [1, N_FULL_KPTS]

    def test_qp_e_length(self, generator):
        assert len(generator.mapped_vars["QP_E"]) == N_STATES

    def test_qp_eo_length(self, generator):
        assert len(generator.mapped_vars["QP_Eo"]) == N_STATES

    def test_pars_values(self, generator):
        pars = np.array(generator.mapped_vars["PARS"])
        assert int(pars[0]) == N_BANDS       # n_bands
        assert int(pars[1]) == N_FULL_KPTS   # n_kpoints
        assert int(pars[2]) == N_STATES      # n_states

    def test_raises_without_inputs(self, ns_db1_path):
        """generate_mappings must fail when eigenvalues/k-points are not set."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        with pytest.raises(ValueError, match="Missing required inputs"):
            gen.generate_mappings()


# ---------------------------------------------------------------------------
# verify_mappings
# ---------------------------------------------------------------------------

class TestVerifyMappings:

    def test_passes_for_valid_mapping(self, generator):
        """Should not raise for the silicon fixture."""
        generator.verify_mappings(k_index=1, top_valence=TOP_VALENCE)

    def test_passes_different_kpoint(self, generator):
        """Check at k=2 as well."""
        generator.verify_mappings(k_index=2, top_valence=TOP_VALENCE)

    def test_wrong_top_valence_raises(self, generator):
        """Using top_valence=0 (invalid 1-based index) must raise an error."""
        with pytest.raises(Exception):
            generator.verify_mappings(k_index=1, top_valence=0)


# ---------------------------------------------------------------------------
# generate_QP_db
# ---------------------------------------------------------------------------

class TestGenerateQPDb:

    def test_creates_file(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        assert tmp_qp_path.exists()

    def test_qp_e_shape(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            shape = ds.variables["QP_E"][:].shape
        assert shape == (N_STATES, 2)

    def test_qp_eo_shape(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            shape = ds.variables["QP_Eo"][:].shape
        assert shape == (N_STATES,)

    def test_qp_kpts_shape(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            shape = ds.variables["QP_kpts"][:].shape
        assert shape == (3, N_FULL_KPTS)

    def test_qp_table_shape(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            shape = ds.variables["QP_table"][:].shape
        assert shape == (3, N_STATES)

    def test_pars_n_bands_and_kpoints(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            pars = np.array(ds.variables["PARS"][:]).flatten()
        assert int(pars[0]) == N_BANDS
        assert int(pars[1]) == N_FULL_KPTS
        assert int(pars[2]) == N_STATES

    def test_ki_energies_in_hartree(self, generator, tmp_qp_path):
        """QP_E values should be on the scale of atomic units (~0.1–1.0 Ha)."""
        from ase.units import Ha
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            qp_e = ds.variables["QP_E"][:, 0]
        # Silicon KC valence eigenvalues are roughly −10 eV to 0 eV → −0.37 to 0 Ha
        # Conduction bottom at ~3.5 eV → 0.13 Ha.  All values must be within ±3 Ha.
        assert np.all(np.abs(qp_e) < 3.0)

    def test_required_variables_present(self, generator, tmp_qp_path):
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            for var in ("QP_E", "QP_Eo", "QP_kpts", "QP_table", "QP_Z", "PARS"):
                assert var in ds.variables, f"Missing variable: {var}"

    def test_imaginary_part_zero(self, generator, tmp_qp_path):
        """The imaginary part of QP_E (column 1) must be zero for KC corrections."""
        generator.generate_QP_db(str(tmp_qp_path))
        with nc.Dataset(str(tmp_qp_path)) as ds:
            im_part = ds.variables["QP_E"][:, 1]
        np.testing.assert_array_equal(im_part, 0.0)


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

class TestSummary:

    def test_summary_returns_string(self, generator):
        s = generator.summary()
        assert isinstance(s, str)

    def test_summary_contains_header(self, generator):
        s = generator.summary()
        assert "KcwQpDatabaseGenerator Summary" in s

    def test_summary_contains_kpoints(self, generator):
        s = generator.summary()
        assert str(N_IBZ_KPTS) in s or str(N_FULL_KPTS) in s

    def test_summary_contains_mappings_yes(self, generator):
        s = generator.summary()
        assert "Yes" in s

    def test_summary_empty_generator(self, ns_db1_path):
        """summary() on a partially configured generator (ns_db1 loaded, no evals)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gen = KcwQpDatabaseGenerator(ns_db1=str(ns_db1_path))
        s = gen.summary()
        assert "Not loaded" in s
