# K2Y Package Improvements

## Overview
This document outlines recommendations for improving the k2y package based on code review.

## 1. Code Quality Improvements

### Current Issues
1. **Inconsistent type hints**: Some methods have them, others don't
2. **Minimal error handling**: No try-except blocks for file I/O operations
3. **Print statements for logging**: Should use Python's `logging` module
4. **Magic numbers**: Hard-coded dimension names like 'D_0000000064'
5. **Long methods**: `generate_mappings()` is >150 lines and does multiple things
6. **Commented-out code**: Several blocks of dead code should be removed

### Recommendations

#### A. Add Proper Logging
```python
import logging

logger = logging.getLogger(__name__)

# Replace print() with logger.info(), logger.warning(), logger.debug()
logger.info("Reading ns.db1 from: %s", ns_db1)
logger.warning("KCW has more bands (%d) than ns.db1 (%d)", n_bands_kcw, n_bands_ns_db1)
```

#### B. Add Type Hints Throughout
```python
def generate_mappings(self) -> None:
    """Generate k-point and eigenvalue mappings."""
    ...

def set_kpoints_from_pwinput(self, pwinput_path: Union[str, Path]) -> Tuple[np.ndarray, str]:
    """Extract k-points from QE input file."""
    ...
```

#### C. Add Error Handling
```python
def __init__(self, ns_db1: Optional[Union[str, Path]] = None, ...):
    if ns_db1:
        try:
            self.ns_db1 = xarray.open_dataset(ns_db1, engine='netcdf4')
        except FileNotFoundError:
            raise FileNotFoundError(f"ns.db1 file not found: {ns_db1}")
        except Exception as e:
            raise RuntimeError(f"Failed to load ns.db1: {e}")
```

#### D. Extract Helper Methods
Break down `generate_mappings()` into smaller, testable functions:

```python
def generate_mappings(self) -> None:
    """Generate k-point and eigenvalue mappings."""
    self._prepare_eigenvalues()
    self._expand_kpoint_grid()
    self._map_kpoints_to_yambo_grid()
    self._compute_qp_corrections()
    self._prepare_output_variables()

def _prepare_eigenvalues(self) -> Tuple[np.ndarray, np.ndarray]:
    """Handle spin dimensions and reshape eigenvalues if needed."""
    ...

def _map_kpoints_to_yambo_grid(self) -> Dict[int, int]:
    """Create mapping between KCW and Yambo k-point grids."""
    ...
```

#### E. Use Configuration/Constants
```python
# At module level or in a config file
QP_DATABASE_DIMENSIONS = {
    'QP_EO_DIM': 'D_0000000064',
    'QP_E_DIM': 'D_0000000077',
    'KPTS_DIM_X': 'D_0000000008',
    'KPTS_DIM_Y': 'D_0000000007',
}

DEFAULT_TOLERANCE = 1e-4  # For k-point matching
ENERGY_UNIT_CONVERSION = Ha  # ASE Hartree to eV
```

## 2. Documentation Improvements

### Add Comprehensive Docstrings

**Current state**: Most methods have minimal or no docstrings  
**Recommendation**: Follow NumPy docstring style

Example:
```python
def generate_mappings(self) -> None:
    """
    Generate k-point and eigenvalue mappings between KCW and Yambo grids.
    
    This method performs the core conversion logic:
    1. Handles spin dimensions in ns.db1 (if present)
    2. Reshapes flattened eigenvalue arrays
    3. Expands k-point grids to full BZ
    4. Maps k-points between KCW and Yambo grids
    5. Computes QP corrections: ΔE = E_KI + (E_KS^KCW - E_KS^Yambo)
    
    The QP correction formula ensures that when Yambo applies:
        E_QP = E_KS^Yambo + ΔE
    it yields the Koopmans eigenvalues E_KI.
    
    Raises
    ------
    ValueError
        If eigenvalue shapes are incompatible
    RuntimeError
        If k-point mapping fails
        
    Notes
    -----
    - Requires `eigenvalues_KI` and `eigenvalues_KS` to be set
    - Requires `kpoints_grid_kcw` to be set
    - Updates attributes: `new_KI`, `new_KS`, `mapped_vars`, `mapped_dims`
    
    See Also
    --------
    set_koopmans_eval : Load Koopmans eigenvalues
    set_kpoints_from_pwinput : Load k-points from QE input
    verify_mappings : Validate the generated mappings
    """
```

### Add README Sections

Update `/home/aiida/data/codes/k2y/README.md`:

```markdown
## Installation

```bash
pip install -e .  # From source
```

## Quick Start

```python
from k2y import KcwQpDatabaseGenerator

# Initialize converter
converter = KcwQpDatabaseGenerator(ns_db1="SAVE/ns.db1")

# Load Koopmans eigenvalues
converter.set_koopmans_eval(path="kc.kho")

# Load k-points
converter.set_kpoints_from_pwinput("nscf.in")

# Generate mappings
converter.generate_mappings()

# Verify (optional but recommended)
converter.verify_mappings(k_index=1, top_valence=10)

# Generate QP database
converter.generate_QP_db("ndb.QP")
```

## API Reference

See [API Documentation](docs/api.md) for detailed reference.

## Troubleshooting

### Common Issues

**Q: "IndexError: index 1 is out of bounds for axis 0 with size 1"**  
A: This occurs when ns.db1 has spin dimensions. The code now handles this automatically.

**Q: "ValueError: operands could not be broadcast together"**  
A: Eigenvalue arrays have incompatible shapes. Ensure KCW and Yambo calculations use compatible k-point meshes.

**Q: K-point mapping fails**  
A: Check that k-points in pwnscf.in match those used in the KCW calculation.
```

## 3. Testing Improvements

### Add Unit Tests

Create `/home/aiida/data/codes/k2y/tests/test_k2y.py`:

```python
import pytest
import numpy as np
from k2y import KcwQpDatabaseGenerator

def test_init():
    """Test initialization without files."""
    converter = KcwQpDatabaseGenerator()
    assert converter.ns_db1 is None
    assert converter.eigenvalues_KI is None

def test_reshape_eigenvalues():
    """Test eigenvalue reshaping."""
    converter = KcwQpDatabaseGenerator()
    
    # Simulate flattened eigenvalues
    n_k, n_b = 64, 200
    flat_evals = np.random.rand(n_k * n_b)
    
    converter.eigenvalues_KS = flat_evals
    converter.kpoints_grid_kcw = np.random.rand(n_k, 3)
    
    # Should reshape to (n_k, n_b)
    # Test this method when extracted

def test_kpoint_matching_tolerance():
    """Test k-point matching with numerical tolerance."""
    k1 = np.array([0.0, 0.0, 0.0])
    k2 = np.array([1e-6, 1e-6, 1e-6])
    
    # Should match within tolerance
    assert np.all(np.abs(k1 - k2) < 1e-4)

@pytest.mark.integration
def test_full_workflow_with_mock_data():
    """Integration test with mock data."""
    # Create mock ns.db1, eigenvalues, etc.
    # Run full workflow
    # Verify output file structure
    pass
```

### Add Integration Tests with Real Data

```python
@pytest.mark.slow
@pytest.mark.requires_data
def test_silicon_conversion():
    """Test conversion with real Silicon data."""
    # Use example data from examples/
    ...
```

## 4. Performance Improvements

### A. Vectorize K-point Matching

**Current**: Loop over all k-points  
**Improvement**: Use broadcasting and vectorized operations

```python
# Current (slow for large grids)
for k in range(n_kpoints_yambo):
    where_first = np.where(np.all(np.abs(self.kpoints_grid_kcw - self.yambopy_ns_db1.red_kpoints[k])<1e-4, axis=1))[0]
    ...

# Improved (vectorized)
def find_kpoint_mapping(kpoints_a, kpoints_b, tolerance=1e-4):
    """Vectorized k-point matching."""
    # Use cdist or similar for efficient distance computation
    from scipy.spatial.distance import cdist
    distances = cdist(kpoints_a, kpoints_b)
    matches = np.argmin(distances, axis=1)
    valid = distances[np.arange(len(kpoints_a)), matches] < tolerance
    return matches, valid
```

### B. Lazy Loading

```python
@property
def ns_db1(self):
    """Lazy load ns.db1 only when needed."""
    if self._ns_db1 is None:
        if self._ns_db1_path is not None:
            self._ns_db1 = xarray.open_dataset(self._ns_db1_path, engine='netcdf4')
    return self._ns_db1
```

## 5. User Experience Improvements

### A. Progress Bars

For long operations:

```python
from tqdm import tqdm

for k in tqdm(range(n_kpoints_yambo), desc="Mapping k-points"):
    ...
```

### B. Validation Methods

```python
def validate_inputs(self) -> Dict[str, bool]:
    """
    Validate that all required inputs are set.
    
    Returns
    -------
    dict
        Dictionary of validation results
        
    Examples
    --------
    >>> converter.validate_inputs()
    {
        'ns_db1': True,
        'eigenvalues_KI': True,
        'eigenvalues_KS': True,
        'kpoints_grid_kcw': False  # Missing!
    }
    """
    return {
        'ns_db1': self.ns_db1 is not None,
        'eigenvalues_KI': self.eigenvalues_KI is not None,
        'eigenvalues_KS': self.eigenvalues_KS is not None,
        'kpoints_grid_kcw': self.kpoints_grid_kcw is not None,
    }

def validate_or_raise(self) -> None:
    """Validate inputs and raise informative error if missing."""
    validation = self.validate_inputs()
    missing = [k for k, v in validation.items() if not v]
    if missing:
        raise ValueError(f"Missing required inputs: {', '.join(missing)}")
```

### C. Summary Method

```python
def summary(self) -> str:
    """
    Generate a summary of the current state.
    
    Returns
    -------
    str
        Multi-line summary string
    """
    lines = [
        "KcwQpDatabaseGenerator Summary",
        "=" * 40,
        f"ns.db1 loaded: {self.ns_db1 is not None}",
    ]
    
    if self.ns_db1 is not None:
        lines.append(f"  K-points in ns.db1: {self.yambopy_ns_db1.nkpoints}")
        lines.append(f"  Bands in ns.db1: {self.ns_db1_evalues.shape[1] if hasattr(self, 'ns_db1_evalues') else 'N/A'}")
    
    if self.eigenvalues_KI is not None:
        lines.append(f"Koopmans eigenvalues: {self.eigenvalues_KI.shape}")
    
    if self.kpoints_grid_kcw is not None:
        lines.append(f"KCW k-points: {self.kpoints_grid_kcw.shape[0]} ({self.kpoints_type})")
    
    if self.mapped_vars is not None:
        lines.append("Mappings generated: Yes")
        lines.append(f"  Output bands: {self.mapped_vars['QP_QP_@_state_1_b_range']}")
        lines.append(f"  Output k-points: {self.mapped_vars['QP_QP_@_state_1_K_range']}")
    
    return "\n".join(lines)
```

## 6. Project Structure

### Suggested Organization

```
k2y/
├── src/
│   └── k2y/
│       ├── __init__.py
│       ├── core.py              # Main KcwQpDatabaseGenerator class
│       ├── utils/
│       │   ├── __init__.py
│       │   ├── kpoints.py       # K-point utilities
│       │   ├── eigenvalues.py   # Eigenvalue handling
│       │   └── io.py            # I/O utilities
│       ├── aiida_interface.py   # AiiDA-specific code
│       ├── cli.py               # Command-line interface
│       └── templates/
├── tests/
│   ├── unit/
│   ├── integration/
│   └── data/                    # Test data
├── examples/
├── docs/
│   ├── api.md
│   ├── tutorial.md
│   └── theory.md                # Explain the physics/math
├── pyproject.toml
└── README.md
```

## 7. Priority Improvements

### High Priority (Do First)
1. ✅ Add comprehensive docstrings to all public methods
2. ✅ Add type hints throughout
3. ✅ Replace print() with logging
4. ✅ Add input validation with clear error messages
5. Add unit tests for critical functions

### Medium Priority
6. Extract helper methods from `generate_mappings()`
7. Add progress indicators for long operations
8. Improve README with troubleshooting section
9. Add `summary()` and `validate_inputs()` methods

### Low Priority
10. Vectorize k-point matching for performance
11. Add integration tests
12. Create comprehensive API documentation
13. Add CLI improvements

## 8. Backward Compatibility

When making improvements, ensure:
- Existing scripts continue to work
- Add deprecation warnings for removed features
- Maintain existing method signatures (add new optional parameters)

Example:
```python
def set_koopmans_eval(
    self, 
    path: Optional[Union[str, Path]] = None,
    output_ase=None,  # Keep for backward compatibility
    results: str = "pki_eigenvalues_on_grid"
) -> None:
    """Load Koopmans eigenvalues."""
    if output_ase is not None:
        warnings.warn(
            "Parameter 'output_ase' is deprecated. Use 'path' instead.",
            DeprecationWarning,
            stacklevel=2
        )
```

## Summary

The k2y package is functional but could benefit from:
1. **Better documentation** (docstrings, README, theory docs)
2. **More robust error handling** (validation, clear error messages)
3. **Better code organization** (extract helpers, use logging)
4. **Testing infrastructure** (unit + integration tests)
5. **Performance improvements** (vectorization for large systems)

These improvements will make the package more maintainable, easier to use, and more reliable for production workflows.
