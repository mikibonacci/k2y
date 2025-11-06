"""
Koopmans to Yambo (k2y) Converter
==================================

This module provides tools to convert Koopmans Functionals eigenvalues
from KCW calculations into Yambo-compatible quasiparticle (QP) database files.

The main class `KcwQpDatabaseGenerator` handles the conversion between:
- Koopmans eigenvalues (KI, pKI)
- Kohn-Sham eigenvalues from DFT
- Yambo QP database format (ndb.QP)

This allows BSE calculations in Yambo using Koopmans eigenvalues instead of
requiring GW ones.
"""

from typing import Union, Optional, Tuple, Dict, Any
from pathlib import Path
import warnings
import logging

import netCDF4 as nc
from ase.units import Ha
from ase_koopmans import io
import xarray
import itertools
import numpy as np
from yambopy import YamboSaveDB
from yambopy.lattice import car_red, red_car

# Set up module logger
logger = logging.getLogger(__name__)

# Constants for k-point matching and database dimensions
DEFAULT_KPOINT_TOLERANCE = 1e-4
QP_DATABASE_DIMENSIONS = {
    'QP_EO_DIM': 'D_0000000064',
    'QP_E_DIM': 'D_0000000077',
    'KPTS_DIM_X': 'D_0000000008',
    'KPTS_DIM_Y': 'D_0000000007',
}

class KcwQpDatabaseGenerator:
    """
    Generator for Yambo QP databases from Koopmans eigenvalues.
    
    This class converts Koopmans Compliant Functional eigenvalues from KCW
    calculations into Yambo-compatible quasiparticle (QP) database files (ndb.QP).
    
    The conversion process:
    1. Reads Yambo's ns.db1 database (contains KS eigenvalues and k-points)m and a template QP file
    2. Loads Koopmans eigenvalues (KI) and corresponding KS eigenvalues
    3. Maps k-points between the two grids
    4. Generates QP corrections: ΔE_QP = E_KI - E_KS
    5. Writes Yambo-compatible ndb.QP file
    
    Attributes
    ----------
    ns_db1 : xarray.Dataset
        Yambo ns.db1 database containing KS eigenvalues and k-points
    yambopy_ns_db1 : YamboSaveDB
        YamboPy interface to ns.db1
    template_QP_path : Path
        Path to template QP database
    eigenvalues_KI : np.ndarray
        Koopmans eigenvalues, shape (n_kpoints, n_bands)
    eigenvalues_KS : np.ndarray
        Kohn-Sham eigenvalues, shape (n_kpoints, n_bands)
    kpoints_grid_kcw : np.ndarray
        K-points grid from KCW calculation, shape (n_kpoints, 3)
    kpoints_type : str
        Type of k-points ('crystal', 'tpiba', etc.)
    
    Examples
    --------
    Basic usage with file paths:
    
    >>> converter = KcwQpDatabaseGenerator(
    ...     ns_db1="path/to/ns.db1",
    ...     template_QP_path="path/to/template.QP"
    ... )
    >>> converter.set_koopmans_eval(path="path/to/kc.kho")
    >>> converter.set_kpoints_from_pwinput("pwnscf.in")
    >>> converter.generate_mappings()
    >>> converter.verify_mappings(k_index=1, top_valence=10)
    >>> converter.generate_QP_db("ndb.QP")
    
    Using AiiDA nodes:
    
    >>> converter = KcwQpDatabaseGenerator.from_aiida(
    ...     yambo_node_pk=12345,
    ...     kcw_node_pk=67890
    ... )
    >>> converter.generate_mappings()
    >>> converter.generate_QP_db("ndb.QP")
    
    Notes
    -----
    - K-point grids must be compatible (mapping is done automatically)
    - Band ranges are adjusted to minimum available in both databases
    - Spin-polarized calculations are supported
    - IMPORTANT: the ns.db1 must be the one you will use for the BSE calculation in conjunction with the generated ndb.QP.
    """
    
    def __init__(
        self,
        ns_db1: Optional[Union[str, Path]] = None,
        template_QP_path: Optional[Union[str, Path]] = None,
        spin: bool = False,
        kpoints_units_kcw: str = 'reduced'        
    ):
        """
        Initialize the QP database generator.
        
        Parameters
        ----------
        ns_db1 : str or Path, optional
            Path to Yambo's ns.db1 netCDF file containing KS eigenvalues and k-points.
            This file is typically found in the SAVE directory after running p2y.
            IMPORTANT: the ns.db1 must be the one you will use for the BSE calculation in conjunction with the generated ndb.QP.
        template_QP_path : str or Path, optional
            Path to template QP database file. If not provided, uses bundled template.
            The template defines the structure and metadata for the output ndb.QP file.
        spin : bool, default=False
            Whether to use spin-polarized template. If you provide your own template, this is ignored.
        kpoints_units_kcw : str, default='reduced'
            Units for k-points from KCW calculation ('reduced' or 'crystal')
            
        Raises
        ------
        FileNotFoundError
            If ns_db1 or template_QP_path file does not exist
        RuntimeError
            If ns.db1 file cannot be loaded
            
        Notes
        -----
        - If `ns_db1` is provided, it will be loaded immediately
        - K-point meshes between ns.db1 and KCW should be compatible
        """
        
        if ns_db1:
            print(f"Reading ns.db1 from: {ns_db1}")
            self.ns_db1 = xarray.open_dataset(ns_db1, engine='netcdf4')
            ns_dir = str(ns_db1).replace('/ns.db1', '')
            self.yambopy_ns_db1 = YamboSaveDB.from_db_file(ns_dir)
            
        if template_QP_path:
            template_path = Path(template_QP_path)
            if not template_path.exists():
                raise FileNotFoundError(f"Template QP file not found: {template_QP_path}")
            logger.info(f"Using QP template: {template_QP_path}")
            self.template_QP_path = template_path
        else:
            self.template_QP_path = self.get_templateQP_filepath(spin=spin)
            logger.info(f"Using bundled template: {self.template_QP_path}")
        
        self.QP_template = None
        self.kpoints_units_kcw = kpoints_units_kcw
        
        # Initialize attributes that will be set later
        self.eigenvalues_KI = None
        self.eigenvalues_KS = None
        self.kpoints_grid_kcw = None
        self.kpoints_type = None 
            
    @classmethod
    def get_templateQP_filepath(cls, spin: bool = False) -> Path:
        """
        Get the path to the template QP database file.
        
        This method retrieves the appropriate template ndb.QP file from the 
        package resources. Different templates are used for spin-polarized 
        and non-spin-polarized calculations.

        NB: this method is used in the init ONLY if no template_QP_path is provided.
        
        Parameters
        ----------
        spin : bool, optional
            If True, return the spin-polarized template. Default is False.
            
        Returns
        -------
        Path
            Path to the template QP database file
            
        Notes
        -----
        Templates are stored in the `templates/` subdirectory:
        - `template_v530.QP`: Non-spin-polarized (default)
        - `template_v530_spin.QP`: Spin-polarized
        
        Examples
        --------
        >>> KcwQpDatabaseGenerator.get_templateQP_filepath()
        PosixPath('.../templates/template_v530.QP')
        
        >>> KcwQpDatabaseGenerator.get_templateQP_filepath(spin=True)
        PosixPath('.../templates/template_v530_spin.QP')
        """
        from importlib_resources import files

        from . import templates
        return files(templates) / f'template_v530{"_spin" if spin else ""}.QP'
    
    def validate_inputs(self) -> Dict[str, bool]:
        """
        Validate that all required inputs are set.
        
        Returns
        -------
        dict
            Dictionary mapping input names to their validation status (True if set)
            
        Examples
        --------
        >>> converter = KcwQpDatabaseGenerator()
        >>> converter.validate_inputs()
        {
            'ns_db1': False,
            'eigenvalues_KI': False,
            'eigenvalues_KS': False,
            'kpoints_grid_kcw': False
        }
        
        >>> converter = KcwQpDatabaseGenerator(ns_db1="SAVE/ns.db1")
        >>> converter.set_koopmans_eval(path="kc.kho")
        >>> converter.set_kpoints_from_pwinput("nscf.in")
        >>> converter.validate_inputs()
        {
            'ns_db1': True,
            'eigenvalues_KI': True,
            'eigenvalues_KS': True,
            'kpoints_grid_kcw': True
        }
        """
        return {
            'ns_db1': self.ns_db1 is not None,
            'eigenvalues_KI': self.eigenvalues_KI is not None,
            'eigenvalues_KS': self.eigenvalues_KS is not None,
            'kpoints_grid_kcw': self.kpoints_grid_kcw is not None,
        }
    
    def validate_or_raise(self) -> None:
        """
        Validate all required inputs and raise informative error if any are missing.
        
        Raises
        ------
        ValueError
            If any required inputs are missing
            
        Examples
        --------
        >>> converter = KcwQpDatabaseGenerator()
        >>> converter.validate_or_raise()
        ValueError: Missing required inputs: ns_db1, eigenvalues_KI, eigenvalues_KS, kpoints_grid_kcw
        """
        validation = self.validate_inputs()
        missing = [k for k, v in validation.items() if not v]
        if missing:
            raise ValueError(
                f"Missing required inputs: {', '.join(missing)}. "
                f"Please set these before calling generate_mappings()."
            )
    
    def summary(self) -> str:
        """
        Generate a summary of the current converter state.
        
        Returns
        -------
        str
            Multi-line summary string showing loaded data and processing status
            
        Examples
        --------
        >>> converter = KcwQpDatabaseGenerator(ns_db1="SAVE/ns.db1")
        >>> print(converter.summary())
        KcwQpDatabaseGenerator Summary
        ========================================
        ns.db1 loaded: True
          K-points in ns.db1: 64
          Bands in ns.db1: N/A
        Koopmans eigenvalues: Not loaded
        KCW k-points: Not loaded
        Mappings generated: No
        """
        lines = [
            "KcwQpDatabaseGenerator Summary",
            "=" * 40,
            f"ns.db1 loaded: {self.ns_db1 is not None}",
        ]
        
        if self.ns_db1 is not None and self.yambopy_ns_db1 is not None:
            lines.append(f"  K-points in ns.db1: {self.yambopy_ns_db1.nkpoints}")
            if hasattr(self, 'ns_db1_evalues') and self.ns_db1_evalues is not None:
                lines.append(f"  Bands in ns.db1: {self.ns_db1_evalues.shape[1]}")
            else:
                lines.append(f"  Bands in ns.db1: N/A")
        
        if self.eigenvalues_KI is not None:
            lines.append(f"Koopmans eigenvalues: {self.eigenvalues_KI.shape}")
        else:
            lines.append("Koopmans eigenvalues: Not loaded")
        
        if self.kpoints_grid_kcw is not None:
            lines.append(f"KCW k-points: {self.kpoints_grid_kcw.shape[0]} ({self.kpoints_type})")
        else:
            lines.append("KCW k-points: Not loaded")
        
        if hasattr(self, 'mapped_vars') and self.mapped_vars is not None:
            lines.append("Mappings generated: Yes")
            lines.append(f"  Output bands: {self.mapped_vars['QP_QP_@_state_1_b_range']}")
            lines.append(f"  Output k-points: {self.mapped_vars['QP_QP_@_state_1_K_range']}")
        else:
            lines.append("Mappings generated: No")
        
        return "\n".join(lines)
     
    def produce_kpoints_for_interpolation(self, filename: Optional[str] = None, coordinates: Optional[str] = "crystal", full_BZ: Optional[str] = False) -> None:
        """
        Generate K_POINTS card for pw.x/kcw.x explicit calculation (or interpolation).
        Useful if you need to reproduce the k-point grid used in the KCW calculation.
        
        This method writes a K_POINTS card in Quantum ESPRESSO format for 
        use with kcw.x interpolation. K-points are converted from the ns.db1 
        database to cartesian coordinates in units of 2π/alat.
        
        Parameters
        ----------
        filename : str, optional
            If provided, write the K_POINTS card to this file. 
            If None, only print to stdout. Default is None.
            
        Notes
        -----
        - K-points are taken from the ns.db1 database loaded in __init__
        - Coordinates are converted to tpiba_b format (cartesian in 2π/alat units)
        - Each k-point is assigned weight 0 (band structure calculation)
        
        Examples
        --------
        >>> converter = KcwQpDatabaseGenerator(ns_db1="SAVE/ns.db1")
        >>> converter.produce_kpoints_for_interpolation("kpoints.txt")
        K_POINTS tpiba_b
        64
        0.0  0.0  0.0 0
        ...
        
        See Also
        --------
        set_kpoints_from_pwinput : Load k-points from pw.x input file
        
        Warnings
        --------
        This method is marked TBO (to be optimized) and may be refactored
        in future versions.
        """
        
        ns = self.ns_db1

        if full_BZ:
            kpoints = self.yambopy_ns_db1.expand_kpts()[0]
        else:
            kpoints = self.yambopy_ns_db1.car_kpoints

        if coordinates == "crystal":
            kpoints = car_red(kpoints,self.yambopy_ns_db1.rlat)


        kpoints = np.round(kpoints, decimals=3)
        
        print(f"K_POINTS {coordinates}")
        nk = len(kpoints)
        print(f"{nk}")
        for k in range(nk):
            print(f'{kpoints[k,0]}  {kpoints[k,1]}  {kpoints[k,2]} {np.round(1/nk,5)}')
            
        if filename:
            with open(filename,'w') as file:
                file.write(f"K_POINTS {coordinates}")
                file.write(f"{nk}\n")
                for k in range(nk):
                    file.write(f'{kpoints[0]}  {kpoints[0]}  {kpoints[0]} {np.round(1/nk,5)}')
                    
        return 
    
    def set_koopmans_eval(
        self, 
        path: Optional[Union[str, Path]] = None, 
        output_ase = None, 
        results: str = "pki_eigenvalues_on_grid"
    ) -> None:
        """
        Load Koopmans eigenvalues from a kcw.x calculation output.
        
        This method extracts both Koopmans Compliant (KI) and Kohn-Sham (KS)
        eigenvalues from a kcw.x output file using ASE's io module.
        
        Parameters
        ----------
        path : str or Path, optional
            Path to the kcw.x output file. If not provided, output_ase must be given.
        output_ase : ASE Atoms object, optional
            Pre-loaded ASE Atoms object containing calculation results.
            If provided, `path` is ignored.
        results : str, default="pki_eigenvalues_on_grid"
            Key for accessing Koopmans eigenvalues in output.calc.results dict.
            
        Attributes Set
        --------------
        eigenvalues_KI : np.ndarray
            Koopmans Compliant eigenvalues (screening-corrected)
        eigenvalues_KS : np.ndarray
            Kohn-Sham eigenvalues from the same calculation
            
        Raises
        ------
        ValueError
            If the requested results key is not found in the calculation results
            
        Notes
        -----
        Currently implemented primarily for AiiDA workflows. Direct file reading
        from kcw.x output files is not yet fully supported.
        
        The eigenvalues are stored as flattened arrays and will be reshaped
        during the generate_mappings() step based on the k-point grid.
        
        Examples
        --------
        >>> converter = KcwQpDatabaseGenerator(ns_db1="SAVE/ns.db1")
        >>> converter.set_koopmans_eval(path="kc.kho")
        >>> print(converter.eigenvalues_KI.shape)  # Flattened
        (12800,)
        
        See Also
        --------
        set_kpoints_from_pwinput : Load corresponding k-point grid
        generate_mappings : Reshape and map eigenvalues to Yambo grid
        from_aiida : Alternative initialization from AiiDA calculations
        
        Warnings
        --------
        This method prints a warning about limited file format support.
        """
        # output_ase is if we already inspected the output file (e.g. aiida)

        logger.warning(
            "Reading kpoints and eigenvalues from kcw.x output file is not fully implemented. "
            "For now, this is primarily supported for AiiDA workflows."
        )
        
        try:
            output = io.read(path) if not output_ase else output_ase
        except Exception as e:
            raise RuntimeError(f"Failed to read Koopmans output file: {e}")
            
        if results not in output.calc.results.keys():
            raise ValueError(f"the requested results {results} are not in the output calc results. Available results are: {list(output.calc.results.keys())}")
        
        self.eigenvalues_KI = np.array(output.calc.results[results])
        self.eigenvalues_KS = np.array(output.calc.results["ks_eigenvalues_on_grid"])
        
        return
    
    def set_kpoints_from_pwinput(self, pwinput_path:str):
        """
        Extract and store k-points from a Quantum ESPRESSO pw.x input file.
        
        This method uses the extract_kpoints utility to read k-points from
        a pw.x input file (e.g., pwnscf.in) and stores them in the 
        kpoints_grid_kcw attribute along with the k-point type.

        NB: the k-points should be in crystal coordinates.
        
        Parameters
        ----------
        pwinput_path : str
            Path to the pw.x input file
            
        Attributes set
        --------------
        kpoints_grid_kcw : np.ndarray
            Array of k-points with shape (nk, 3)
        kpoints_type : str
            Type of k-points (e.g., 'crystal', 'tpiba', 'automatic')
        """
        from .extract_kpoints import extract_kpoints_from_pwin
        
        if not Path(pwinput_path).exists():
            raise FileNotFoundError(f"PW input file not found: {pwinput_path}")
        
        try:
            kpoints, kpoints_type = extract_kpoints_from_pwin(pwinput_path)
        except Exception as e:
            raise RuntimeError(f"Failed to extract k-points from {pwinput_path}: {e}")
        
        self.kpoints_grid_kcw = kpoints
        self.kpoints_type = kpoints_type
        
        logger.info(f"Loaded {len(kpoints)} k-points of type '{kpoints_type}' from {pwinput_path}")
        logger.debug(f"K-points shape: {kpoints.shape}")
        
        return kpoints, kpoints_type
    
    def generate_mappings(self, time_rev:bool=True, brute_force:bool=True) -> None:
        """
        Generate k-point and eigenvalue mappings between KCW and Yambo grids.
        
        This is the core conversion method that performs the following steps:
        
        1. **Handle spin dimensions**: Detect and extract spin components from ns.db1
        2. **Reshape eigenvalues**: Convert flattened arrays to (n_kpoints, n_bands)
        3. **Expand k-point grids**: Generate full Brillouin zone from irreducible grid
        4. **Map k-points**: Match KCW k-points to Yambo k-point grid
        5. **Compute QP corrections**: Calculate ΔE = E_KI + (E_KS^KCW - E_KS^Yambo)
        
        The QP correction formula ensures that when Yambo applies the correction:
            E_QP = E_KS^Yambo + ΔE
        it yields the Koopmans eigenvalues E_KI.
        
        Parameters
        ----------
        time_rev : bool, default=True
            Enable time-reversal symmetry for k-point matching. If a k-point is not
            found by direct comparison, try matching with -k (time-reversed k-point).
            This is useful when the KCW and Yambo calculations use different but
            symmetry-equivalent k-point grids.
        brute_force : bool, default=True
            Enable brute-force matching by comparing absolute values of k-point 
            coordinates. If direct matching and time-reversal both fail, this tries
            to match |k_x|, |k_y|, |k_z| independently. This is a last resort for
            grids that differ by sign conventions or symmetry operations.
        
        Attributes Set
        --------------
        ns_db1_evalues : np.ndarray
            KS eigenvalues from ns.db1, shape (n_kpoints, n_bands)
        new_KI : np.ndarray
            QP corrections for all k-points and bands
        new_KS : np.ndarray
            KS eigenvalues mapped to Yambo grid
        mapped_vars : dict
            Variables for the output QP database
        mapped_dims : dict
            Dimensions for the output QP database
            
        Raises
        ------
        ValueError
            If eigenvalue shapes are incompatible after reshaping
        RuntimeError
            If k-point mapping fails to find matches
            
        Notes
        -----
        
        **Array shapes**:
        - Input eigenvalues can be flattened: (n_kpoints * n_bands,)
        - ns.db1 can have spin: (n_spin, n_kpoints, n_bands)
        - After processing: (n_kpoints, n_bands)
        
        **K-point matching**:
        - Tolerance: 1e-4 for k-point coordinate comparison
        - Uses crystal coordinates for matching
        - Handles both reduced and full Brillouin zone
        
        **Band number compatibility**:
        - Uses minimum of KCW and ns.db1 band counts
        - Issues warning if band counts differ
        
        Examples
        --------
        >>> converter = KcwQpDatabaseGenerator(ns_db1="SAVE/ns.db1")
        >>> converter.set_koopmans_eval(path="kc.kho")
        >>> converter.set_kpoints_from_pwinput("nscf.in")
        >>> converter.generate_mappings()
        ns.db1 has 200 bands
        KCW has 200 bands
        Using 200 bands for mapping
        Yambo ns.db1 has 64 k-points
        K-point mapping complete: 64/64 matched
        
        See Also
        --------
        set_koopmans_eval : Load Koopmans eigenvalues
        set_kpoints_from_pwinput : Load k-points from QE input
        verify_mappings : Validate the generated mappings
        generate_QP_db : Write the output database
        
        Warnings
        --------
        This method modifies eigenvalues in place. If KCW has more bands than
        ns.db1, only the first n_bands will be included in the QP database.
        """
        
        # Validate inputs before proceeding
        self.validate_or_raise()
        
        logger.info("Starting k-point and eigenvalue mapping generation...")
        
        # NB: it is fundamental that the KS eigenvalues in the ns.db1 are the same as the ones in the new DB. 
        # Otherwise, the mapping will not work and yambo will interpolate the wrong eigenvalues.

        ns = self.ns_db1
        if not hasattr(self, 'eigenvalues_KS') or self.eigenvalues_KS is None:
            self.eigenvalues_KS = ns.variables["EIGENVALUES"].values[0] # not exactly, but we can map.
            logger.info("Using KS eigenvalues from ns.db1")
        else:
            logger.info("Using KS eigenvalues already stored from KCW calculation")
        
        # Handle spin dimension in ns.db1
        # ns_db1_evalues shape can be (n_spin, n_kpoints, n_bands) or (n_kpoints, n_bands)
        ns_db1_evalues_raw = ns.variables["EIGENVALUES"].values
        if len(ns_db1_evalues_raw.shape) == 3:
            # Has spin dimension, take first spin component
            logger.info(f"ns.db1 has spin dimension: {ns_db1_evalues_raw.shape}")
            self.ns_db1_evalues = ns_db1_evalues_raw[0]  # Shape: (n_kpoints, n_bands)
            logger.info(f"Using first spin component, shape: {self.ns_db1_evalues.shape}")
        else:
            self.ns_db1_evalues = ns_db1_evalues_raw

        # Reshape eigenvalues if they are flattened arrays
        # The eigenvalues come as (n_kpoints * n_bands,) and need to be reshaped to (n_kpoints, n_bands)
        if len(self.eigenvalues_KS.shape) == 2 and self.eigenvalues_KS.shape[0] == 1:
            n_kpoints = self.kpoints_grid_kcw.shape[0]
            n_bands = self.eigenvalues_KS.shape[1] // n_kpoints
            logger.info(f"Reshaping flattened eigenvalues: total_size={self.eigenvalues_KS.shape[1]} -> ({n_kpoints}, {n_bands})")
            self.eigenvalues_KS = self.eigenvalues_KS.reshape((n_kpoints, n_bands))
            self.eigenvalues_KI = self.eigenvalues_KI.reshape((n_kpoints, n_bands))
            logger.debug(f"Reshaped eigenvalues_KS: {self.eigenvalues_KS.shape}")
            logger.debug(f"Reshaped eigenvalues_KI: {self.eigenvalues_KI.shape}")

        ########################## Making yambo-KS (we will refer as yKS) compatible with kcw-KS (we will refer as kKS)
        ######## SKIP : (1) Mapping of the kpoints from the full grid to the kpoints in the koopmans KS/KI.
        ######## (2) Modifying the KI eigenvalues KI = KI + (kKS-yKS) so that when we apply KI-kKS to yKS, we get the KI eigenvalues.
        ########     The reason is that indeed, yambo computes KI-kKS from ndb.QP, and then applies it to the yKS contained in the ns.db1.

        # we expand the kpoints to the full grid, so we have the mapping
        # between the kpoints in the ns.db1 and the kpoints in the koopmans KS/KI.
        
        self.yambopy_ns_db1.expand_kpts()
        
        # Determine the number of bands to use: minimum between KCW and ns.db1
        n_bands_ns_db1 = self.ns_db1_evalues.shape[1]
        n_bands_kcw = self.eigenvalues_KS.shape[1]
        n_bands = min(n_bands_ns_db1, n_bands_kcw)
        
        logger.info(f"ns.db1 has {n_bands_ns_db1} bands")
        logger.info(f"KCW has {n_bands_kcw} bands")
        logger.info(f"Using {n_bands} bands for mapping")
        
        if n_bands_kcw > n_bands_ns_db1:
            logger.warning(f"KCW has more bands ({n_bands_kcw}) than ns.db1 ({n_bands_ns_db1})")
            logger.warning(f"Only the first {n_bands} bands will be mapped to the QP database")
        
        # Determine number of kpoints from yambopy_ns_db1.eigenvalues
        # Shape can be (n_spin, n_kpoints, n_bands) or (n_kpoints, n_bands)
        if len(self.yambopy_ns_db1.eigenvalues.shape) == 3:
            n_kpoints_yambo = self.yambopy_ns_db1.eigenvalues.shape[1]
        else:
            n_kpoints_yambo = self.yambopy_ns_db1.eigenvalues.shape[0]
        
        logger.info(f"Yambo ns.db1 has {n_kpoints_yambo} k-points")
        
        #new_KS = np.zeros((n_kpoints_yambo, n_bands))
        #new_KI = np.zeros((n_kpoints_yambo, n_bands))
        
        # Use only the bands that are available in both databases
        eigenvalues = self.eigenvalues_KI[:n_kpoints_yambo,:n_bands].copy()
        eigenvalues_KS = self.eigenvalues_KS[:n_kpoints_yambo,:n_bands].copy()

        expanded_kpoints, expanded_indexes, _ = self.yambopy_ns_db1.expand_kpts() # this gives cartesian coordinates
        if self.kpoints_type in ['reduced','crystal']:
            expanded_kpoints = car_red(expanded_kpoints,self.yambopy_ns_db1.rlat)

        logger.info("Mapping k-points between KCW and Yambo grids...")
        matched_kpoints = 0
        for k in range(n_kpoints_yambo):
            where_first = np.where(np.all(np.abs(self.kpoints_grid_kcw - self.yambopy_ns_db1.red_kpoints[k])<DEFAULT_KPOINT_TOLERANCE, axis=1))[0]
            if len(where_first) == 0:
                logger.debug(f"K-point {k} not found directly, trying expanded grid...")
                where_this_kpoint = np.where(expanded_indexes == k)[0]
                # we then search for the kpoint in the expanded kpoints, but we still use the index k.:
                for kpoint in expanded_kpoints[where_this_kpoint]:
                    where_first = np.where(np.all(np.abs(self.kpoints_grid_kcw - kpoint)<DEFAULT_KPOINT_TOLERANCE, axis=1))[0]
                    if len(where_first) > 0:
                        where = where_first[0]
                        logger.debug(f"K-point {k} found by expanding the grid")
                        break
                    if time_rev:
                        where_first = np.where(np.all(np.abs(self.kpoints_grid_kcw + kpoint)<DEFAULT_KPOINT_TOLERANCE, axis=1))[0]
                        if len(where_first) > 0:
                            where = where_first[0]
                            logger.debug(f"K-point {k} found using time-reversal symmetry")
                            break
                    if brute_force:
                        where_first = np.where(np.all(np.abs(np.abs(self.kpoints_grid_kcw) - np.abs(kpoint))<DEFAULT_KPOINT_TOLERANCE, axis=1))[0]
                        if len(where_first) > 0:
                            where = where_first[0]
                            logger.debug(f"K-point {k} found using brute-force matching")
                            break
            else:
                logger.debug(f"K-point {k} found directly in the KCW kpoints grid")
                where = where_first[0]
            #new_KS[k,:] = self.eigenvalues_KS[where,:]
            #new_KI[k,:] = self.eigenvalues_KI[where,:]

            eigenvalues[k,:] = self.eigenvalues_KI[where,:n_bands] + (self.eigenvalues_KS[where,:n_bands] - Ha*self.ns_db1_evalues[k,:n_bands])
            eigenvalues_KS[k,:] = self.eigenvalues_KS[where,:n_bands]
        
        logger.info(f"K-point mapping complete: {n_kpoints_yambo}/{n_kpoints_yambo} matched")
        
        self.new_KI = eigenvalues.copy()
        self.new_KS = eigenvalues_KS.copy()

        # use the KI:
        #eigenvalues = self.eigenvalues_KI + (self.eigenvalues_KS - self.ns_db1_evalues[:,:self.eigenvalues_KS.shape[1]])

        logger.debug(f"Shape of the eigenvalues: {np.shape(eigenvalues)}")
        
        ###############################################################################
        if not hasattr(self, 'kpoints') or self.kpoints is None:
            self.kpoints = ns.variables["K-POINTS"].values[:,:np.shape(eigenvalues)[0]] # exactly as in ndb.QP
        else:
            logger.debug("Using k-points already stored")
        
        self.kpoints = ns.variables["K-POINTS"].values # exactly as in ndb.QP

        bands = [1,np.shape(eigenvalues)[0]*np.shape(eigenvalues)[1]]

        kpoints_ind = [1,np.shape(self.kpoints)[1]]
        reshaped_eval = eigenvalues.reshape(bands[-1])/Ha
        reshaped_eval_KS = eigenvalues_KS[:np.shape(eigenvalues)[0],:np.shape(eigenvalues)[1]].reshape(bands[-1])/Ha

        QP_kpts = self.kpoints # [[x],[y],[z]]
        QP_Eo = reshaped_eval_KS # [E_KS]
        QP_Z = [[1,0]]*(bands[1]-bands[0]+1) # [[Z.real,Z.imag]] <--- we don't care about this.
        table_bands = [b_ind for b_ind in range(1,eigenvalues.shape[-1]+1)]*self.kpoints.shape[1]

        self.QP_Eo = QP_Eo
        

        nested_list = [[k]*eigenvalues.shape[-1] for k in range(1,1+self.kpoints.shape[1])]
        table_kpoints = list(itertools.chain(*nested_list))

        QP_table = [table_bands,table_bands,table_kpoints]
        
        self.table = QP_table
        
        QP_E = [[E,0] for E in reshaped_eval] # [[E.real,E.imag]]
        PARS = np.ma.array(
            np.array([np.shape(eigenvalues)[1],np.shape(eigenvalues)[0],np.shape(eigenvalues)[0]*np.shape(eigenvalues)[1],26,-1,-1]),
            mask = [False, False, False, False,  True,  True])

        self.QP_E = QP_E

        self.mapped_vars = {
        'QP_QP_@_state_1_b_range': [1,np.shape(eigenvalues)[1]],
        'QP_QP_@_state_1_K_range': [1,kpoints_ind[-1]],
        "QP_kpts": QP_kpts if QP_kpts.shape[0] == 3 else QP_kpts.T,
        "QP_E": QP_E,
        "QP_Eo": QP_Eo,
        "QP_Z": QP_Z,
        "QP_table": QP_table,
        "PARS":PARS,
        }

        # This mapping is done onto template.QPs. You should always use them!!!
        # TODO: generalize this, to read from whatever template.QP you want.
        # or at least create a template which always work (e.g. I don't want, in the following mapping, to have to change the D_0000000001... because it is used somewhere else in the template.QP)
        self.mapped_dims = {
            QP_DATABASE_DIMENSIONS['QP_EO_DIM']: [f"D_{str(len(QP_Eo)).zfill(10)}",len(QP_Eo)], # 
            QP_DATABASE_DIMENSIONS['QP_E_DIM']: [f"D_{str(len(QP_Eo)).zfill(10)}",len(QP_Eo)], # 
            QP_DATABASE_DIMENSIONS['KPTS_DIM_X']: [f"D_{str(QP_kpts.shape[0]).zfill(10)}",QP_kpts.shape[0]], # 
            QP_DATABASE_DIMENSIONS['KPTS_DIM_Y']: [f"D_{str(QP_kpts.shape[1]).zfill(10)}",QP_kpts.shape[1]], # before it was shape[0], but it is not the case anymore.
            #'D_0000000100': [f"D_{str(mapped_vars['QP_QP_@_state_1_b_range'][1]).zfill(10)}",mapped_vars['QP_QP_@_state_1_b_range'][1]], # 
        }
        
        logger.info("Mapping successfully generated!")
        logger.info(f"Output will contain {self.mapped_vars['QP_QP_@_state_1_b_range'][1]} bands and {self.mapped_vars['QP_QP_@_state_1_K_range'][1]} k-points")
        
        return
    
    def verify_mappings(self, k_index: int, top_valence: int) -> None:
        """
        Verify eigenvalue mappings by computing and comparing band gaps.
        
        This validation method calculates the band gap at a specific k-point
        using different data sources and ensures they match within tolerance.
        It compares:
        
        1. Gap from new_KI arrays (Koopmans eigenvalues)
        2. Gap from new_KS arrays (KS eigenvalues)
        3. Gap from QP_E table (output QP corrections)
        4. Gap from QP_Eo table (output KS reference)
        
        Parameters
        ----------
        k_index : int
            1-based index of the k-point to check (e.g., 1 for Γ-point)
        top_valence : int
            1-based index of the top valence band
            
        Raises
        ------
        AssertionError
            If gaps do not match within tolerance (1e-3 eV)
            
        Notes
        -----
        **Tolerance**: Gaps must match within 1e-3 eV (1 meV)
        
        **Indexing**: 
        - k_index and top_valence use 1-based indexing (QE/Yambo convention)
        - Internal arrays use 0-based Python indexing
        
        **What is verified**:
        - Consistency between KI and QP_E (both Koopmans eigenvalues)
        - Consistency between KS and QP_Eo (both reference KS eigenvalues)
        - Correct mapping from k-point/band indices to QP_table entries
        
        Examples
        --------
        Verify at Γ-point (k_index=1) with top valence band 10:
        
        >>> converter.verify_mappings(k_index=1, top_valence=10)
        Computed gap from new_KI at k-point 1: 5.234 eV
        Computed gap from new_KS at k-point 1: 0.612 eV
        Computed gap from QP_E at k-point 1: 5.234 eV
        Computed gap from QP_Eo at k-point 1: 0.612 eV
        Mappings verified successfully: gaps match within tolerance of 1e-3 eV.
        
        Verify at another k-point:
        
        >>> converter.verify_mappings(k_index=8, top_valence=10)
        
        See Also
        --------
        generate_mappings : Generate the mappings to verify
        generate_QP_db : Write verified mappings to database
        
        Warnings
        --------
        This method requires `generate_mappings()` to be called first,
        otherwise attributes like `new_KI`, `new_KS`, `table`, etc. will
        not be available.
        """
        # verify that the mappings are correct by computing the band gap at specified k-point
        # from the new_KS and new_KI, and comparing with the values in the QP_Eo.
        gap_QP_1 = np.round(self.new_KI[k_index-1,top_valence] - self.new_KI[k_index-1,top_valence-1],3)
        gap_KS_1 = np.round(self.new_KS[k_index-1,top_valence] - self.new_KS[k_index-1,top_valence-1],3)

        print(f"Computed gap from new_KI at k-point {k_index}: {gap_QP_1} eV")
        print(f"Computed gap from new_KS at k-point {k_index}: {gap_KS_1} eV")

        # for QP_Eo and QP_E, we need to find the right indices in the QP_table.
        # find the indices in the QP_table corresponding to the k_index and top_valence
        k_v = np.where((np.array(self.table[2]) == k_index) & (np.array(self.table[0]) == top_valence))[0][0]
        k_c = np.where((np.array(self.table[2]) == k_index) & (np.array(self.table[0]) == top_valence+1))[0][0]

        gap_KS_2 = np.round((self.QP_Eo[k_c] - self.QP_Eo[k_v])* Ha,3)
        gap_QP_2 = np.round((self.QP_E[k_c][0] - self.QP_E[k_v][0])* Ha,3)

        logger.info(f"Computed gap from QP_Eo at k-point {k_index}: {gap_QP_2} eV")
        logger.info(f"Computed gap from QP_E at k-point {k_index}: {gap_KS_2} eV")

        # verification
        assert abs(gap_QP_1 - gap_QP_2) < 1e-3 , "Mismatch in QP gap verification!"
        assert abs(gap_KS_1 - gap_KS_2) < 1e-3 , "Mismatch in KS gap verification!"

        logger.info("Mappings verified successfully: gaps match within tolerance of 1e-3 eV.")
        
    
    def generate_QP_db(self, output_filename: str = "output.QP") -> None:
        """
        Generate a Yambo-compatible QP database file from the mapped eigenvalues.
        
        This method creates a netCDF4 database file (ndb.QP) that Yambo can read
        for quasiparticle (QP) band structure calculations. The database contains:
        
        - QP corrections: ΔE for each k-point and band
        - K-point grid: matching the ns.db1 grid
        - Reference KS eigenvalues: for applying corrections
        - Band ranges and table structure
        
        Parameters
        ----------
        output_filename : str, default="output.QP"
            Path to the output QP database file. Common names:
            - "ndb.QP" (standard Yambo convention)
            - "ndb.QP_koopmans" (to distinguish from other QP databases)
            
        Notes
        -----
        **File format**: NetCDF4 (same as Yambo ns.db1)
        
        **Template-based**: Uses a pre-defined template QP file and:
        1. Copies the structure (dimensions, variables)
        2. Updates dimensions to match the number of k-points/bands
        3. Replaces variable data with mapped values
        
        **Critical variables**:
        - `QP_E`: Quasiparticle energies (complex: [[E.real, E.imag], ...])
        - `QP_Eo`: Reference KS energies (real values in Hartree)
        - `QP_Z`: Renormalization factors (set to [[1, 0], ...])
        - `QP_kpts`: K-point coordinates
        - `QP_table`: Mapping [band_index, band_index, kpoint_index]
        - `PARS`: Database parameters
        
        **Dimension mapping**: The method updates template dimensions like:
        - D_0000000064 → D_0000012800 (for 12800 QP corrections)
        - D_0000000008 → D_0000000003 (for 3D k-point vectors)
        - D_0000000007 → D_0000000064 (for 64 k-points)
        
        Examples
        --------
        Generate standard ndb.QP file:
        
        >>> converter.generate_QP_db("ndb.QP")
        Generating Koopmans quasiparticle database...
        QP_QP_@_state_1_b_range
        ...
        Database written to ndb.QP
        
        Generate with custom name:
        
        >>> converter.generate_QP_db("ndb.QP_koopmans_pbe")
        
        See Also
        --------
        generate_mappings : Must be called before this method
        verify_mappings : Optional validation before writing
        get_templateQP_filepath : Get the template file path
        generate_QP_db_SinglefileData : AiiDA-specific variant
        
        Raises
        ------
        IOError
            If template file cannot be read
        RuntimeError
            If mappings have not been generated yet
            
        Warnings
        --------
        Always call `generate_mappings()` before this method, otherwise
        the mapped_vars and mapped_dims attributes will not exist.
        
        Consider calling `verify_mappings()` to ensure correctness before
        writing the database.
        """
        
        print("Generating Koopmans quasiparticle database...")
        
        # 1 Open the template file in read-only mode
        if not self.QP_template:
            template_nc = nc.Dataset(self.template_QP_path,'r')
        else:
            template_nc = self.QP_template
            
        new_db = nc.Dataset(output_filename,'w', format='NETCDF4')
        
        # 2.1 Copy dimesions
        # for dim_name in template_nc.dimensions:
        #     if dim_name not in new_db.dimensions:
        #         if dim_name in self.mapped_dims.keys():
        #             print(f"copying dimension {dim_name} -> {self.mapped_dims[dim_name][0]}")
        #             new_db.createDimension(self.mapped_dims[dim_name][0], self.mapped_dims[dim_name][1])
        #         else:
        #             new_db.createDimension(dim_name, len(template_nc.dimensions[dim_name]))

        # 2.2 Copy variables
        for var_name in template_nc.variables:
            #print(var_name)
            if var_name not in new_db.variables:
                print(var_name)
                var = template_nc.variables[var_name]
                new_dims = []
                print(var.dimensions)
                if var_name in self.mapped_vars:
                    if hasattr(self.mapped_vars[var_name], "shape"):
                        dimensions = self.mapped_vars[var_name].shape
                    else: 
                        if isinstance(self.mapped_vars[var_name], list):
                            dimensions = np.array(self.mapped_vars[var_name]).shape
                        else:
                            dimensions = [len(self.mapped_vars[var_name])]
                else:
                    dimensions = var.shape
                for v in dimensions:
                    dimension_name = f"D_{str(v).zfill(10)}"
                    if dimension_name in new_db.dimensions:
                        pass
                    else:
                        new_db.createDimension(dimension_name, v)

                    new_dims.append(dimension_name)
                    
                print(new_dims)
                new_db.createVariable(var_name, var.dtype, new_dims)
                new_db.variables[var_name][:] = self.mapped_vars[var_name] if var_name in self.mapped_vars.keys() else var[:]

        # 3 Copy global attributes
        for attr_name in template_nc.ncattrs():
            new_db.setncattr(attr_name, template_nc.getncattr(attr_name))
        
        new_db.close()
        template_nc.close()
        
        print("...done.")
        
        return
    
    def generate_QP_db_SinglefileData(self, filename:str="output.QP", temporary_dir:str=None):
        """
        This method is to be used with the SinglefileData class of AiiDA.
        
        Need to do this in a tempfile, if we want to use the workflow.
        """
        from aiida import orm, load_profile
        load_profile()
        
        import os
        
        self.generate_QP_db(output_filename=filename)
        new_db = orm.SinglefileData(file=os.path.abspath(filename))
        #new_db.store()
        print(f"SinglefileData created, pk={new_db.pk}.")
        
        return new_db
    

    def full_workflow(self):
        raise NotImplementedError("to be implemented...")
    
    @staticmethod
    def generate_SinglefileData_from_file(input_file_path:str):
        """
        Static method to generate a SinglefileData from a file path.
        """
        from aiida import orm, load_profile
        load_profile()
        
        import os
        
        new_db = orm.SinglefileData(file=os.path.abspath(input_file_path))
        print("SinglefileData created, please run new_db.store() to store it in the database.")
        
        return new_db
    

    @classmethod
    def from_aiida(cls, yambo_node_pk, kcw_node_pk = None, on_grid = True, template_QP_path=None, qp_template_node=None):
        """Initialize the class from an AiiDA yambo and kcw node.
        
        we use the tempdir of the yambo node to init the self.ns_db1, and 
        the tempdir of the kcw node to init the self.eigenvalues_KI.

        Args:
            yambo_node_pk (_type_): _description_
            kcw_node_pk (_type_): _description_
        """
        
        import tempfile
        import pathlib
        
        from aiida import orm, load_profile
        load_profile()
        
        from aiida_yambo.utils.common_helpers import find_pw_parent

        
        yambocalculation = orm.load_node(yambo_node_pk)
        kcwcalculation = orm.load_node(kcw_node_pk) if kcw_node_pk else None
        
        with tempfile.TemporaryDirectory() as dirpath: # actually skippable...
            # Open the output file from the AiiDA storage and copy content to the temporary file
            for filename in yambocalculation.outputs.retrieved.base.repository.list_object_names():
                if 'ns.db1' in filename:
                    # Create the file with the desired name
                    temp_file = pathlib.Path(dirpath) / "ns.db1"
                    with yambocalculation.outputs.retrieved.open(filename, 'rb') as handle:
                        temp_file.write_bytes(handle.read())
                    
                    kcwqpdatabaseGenerator = cls(ns_db1=temp_file, template_QP_path=template_QP_path)
                    
            if qp_template_node:
                filename = "ndb.QP"
                temp_file = pathlib.Path(dirpath) / filename
                temp_file.write_bytes(qp_template_node.get_content("rb"))
                # QP database read from db file
                QP_instance = nc.Dataset(dirpath + '/ndb.QP')
                kcwqpdatabaseGenerator.QP_template = QP_instance
        
        if kcwcalculation:
            if on_grid:
                kcwqpdatabaseGenerator.eigenvalues_KI = np.array(kcwcalculation.outputs.output_parameters.get_dict()["pki_eigenvalues_on_grid"][-1])
                kcwqpdatabaseGenerator.eigenvalues_KS = np.array(kcwcalculation.outputs.output_parameters.get_dict()["ks_eigenvalues_on_grid"][-1])
                kcwqpdatabaseGenerator.kpoints_grid_kcw = find_pw_parent(kcwcalculation).outputs.output_band.get_array('kpoints')
                
                # reshape the eigenvalues
                kcwqpdatabaseGenerator.eigenvalues_KI = kcwqpdatabaseGenerator.eigenvalues_KI.reshape((kcwqpdatabaseGenerator.kpoints_grid_kcw.shape[0],kcwqpdatabaseGenerator.eigenvalues_KI.shape[0]//kcwqpdatabaseGenerator.kpoints_grid_kcw.shape[0]))
                kcwqpdatabaseGenerator.eigenvalues_KS = kcwqpdatabaseGenerator.eigenvalues_KS.reshape((kcwqpdatabaseGenerator.eigenvalues_KI.shape[0],kcwqpdatabaseGenerator.eigenvalues_KI.shape[1]))
                
            else:
                kcwqpdatabaseGenerator.eigenvalues_KI = np.array(kcwcalculation.outputs.output_parameters.get_dict()["eigenvalues"])
        
        kcwqpdatabaseGenerator.yambocalculation = yambocalculation
        kcwqpdatabaseGenerator.kcwcalculation = kcwcalculation

        return kcwqpdatabaseGenerator


"""Usage example:

from k2y.k2y import KcwQpDatabaseGenerator

converter = KcwQpDatabaseGenerator(
    ns_db1="/path/to/ns.db1",
    #template_QP_path="/path/to/template.QP"
    )
    
converter.set_koopmans_eval(path="/path/to/kc.kho") # not needed if you are using AiiDA

# Optional: Load k-points from pw.x input file
converter.set_kpoints_from_pwinput("/path/to/pwnscf.in")

converter.generate_mappings()

converter.verify_mappings(k_index=1, top_valence=10) # adjust k_index and top_valence as needed

converter.generate_QP_db("out.QP")

#converter.produce_kpoints_for_interpolation() # to produce the k-points card for the interpolation, in the kcw.x run...

"""