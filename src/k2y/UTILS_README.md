# K2Y Utilities

This directory contains utility scripts for working with Quantum ESPRESSO and Yambo files.

## extract_kpoints.py

Extract k-points from Quantum ESPRESSO pw.x input files (nscf calculations).

### Features
- Extracts only the first 3 columns (kx, ky, kz coordinates)
- Supports multiple output formats: text, numpy array, Python list
- Handles crystal, tpiba, and other k-point types
- Can output to file or stdout

### Usage

**Print to stdout:**
```bash
python extract_kpoints.py pwnscf.in
```

**Save to text file:**
```bash
python extract_kpoints.py pwnscf.in -o kpoints.txt
```

**Save as numpy array:**
```bash
python extract_kpoints.py pwnscf.in -o kpoints.npy -f numpy
```

**Save as Python list:**
```bash
python extract_kpoints.py pwnscf.in -o kpoints.py -f python
```

**Verbose output:**
```bash
python extract_kpoints.py pwnscf.in -o kpoints.txt -v
```

### Use in Python scripts

You can also import and use the functions directly:

```python
from k2y.extract_kpoints import extract_kpoints_from_pwin, save_kpoints

# Extract k-points
kpoints, kpoints_type = extract_kpoints_from_pwin('pwnscf.in')
print(f"Found {len(kpoints)} k-points of type {kpoints_type}")

# Save in desired format
save_kpoints(kpoints, 'kpoints.txt', format='txt')
```

### Use with KcwQpDatabaseGenerator

The `extract_kpoints` utility is integrated into the `KcwQpDatabaseGenerator` class:

```python
from k2y.k2y import KcwQpDatabaseGenerator

converter = KcwQpDatabaseGenerator(ns_db1="/path/to/ns.db1")

# Load k-points directly from pw.x input file
converter.set_kpoints_from_pwinput("/path/to/pwnscf.in")

# K-points are now stored in:
#   converter.kpoints_grid_kcw  (numpy array of shape (nk, 3))
#   converter.kpoints_type       (string: 'crystal', 'tpiba', etc.)

# Continue with the rest of the workflow
converter.generate_mappings()
converter.generate_QP_db("out.QP")
```

### Output formats

1. **txt** (default): Space-separated text file with 3 columns
   ```
   0.00000000  0.00000000  0.00000000
   0.00000000  0.00000000  0.25000000
   ...
   ```

2. **numpy**: Binary numpy array file (.npy)
   - Can be loaded with `np.load('kpoints.npy')`

3. **python**: Python source file with a list
   ```python
   kpoints = [
       [0.00000000, 0.00000000, 0.00000000],
       [0.00000000, 0.00000000, 0.25000000],
       ...
   ]
   ```
