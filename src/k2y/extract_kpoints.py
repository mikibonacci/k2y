#!/usr/bin/env python3
"""
Utility to extract k-points from Quantum ESPRESSO pw.x nscf input files.

This script reads a pw.x input file and extracts the k-points (first 3 columns only),
which represent the k-point coordinates in crystal or cartesian units.
"""

import argparse
import sys
from pathlib import Path
import numpy as np


def extract_kpoints_from_pwin(filepath, output_format='txt'):
    """
    Extract k-points from a pw.x input file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the pw.x input file
    output_format : str
        Output format: 'txt' (plain text), 'numpy' (npy file), or 'list' (Python list)
    
    Returns
    -------
    kpoints : np.ndarray
        Array of shape (nk, 3) containing k-point coordinates
    kpoints_type : str
        Type of k-points ('crystal', 'tpiba', 'automatic', etc.)
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find K_POINTS card
    kpoints_line_idx = None
    kpoints_type = None
    
    for i, line in enumerate(lines):
        if line.strip().upper().startswith('K_POINTS'):
            kpoints_line_idx = i
            # Extract k-points type (crystal, tpiba, etc.)
            parts = line.strip().split()
            if len(parts) > 1:
                kpoints_type = parts[1].lower()
            else:
                kpoints_type = 'tpiba'  # default
            break
    
    if kpoints_line_idx is None:
        raise ValueError("K_POINTS card not found in input file")
    
    # Next line should contain number of k-points
    nk_line = lines[kpoints_line_idx + 1].strip()
    try:
        nk = int(nk_line.split()[0])
    except (ValueError, IndexError):
        raise ValueError(f"Could not parse number of k-points from line: {nk_line}")
    
    # Extract k-points (first 3 columns only)
    kpoints = []
    for i in range(nk):
        line_idx = kpoints_line_idx + 2 + i
        if line_idx >= len(lines):
            raise ValueError(f"Expected {nk} k-points but file ended early")
        
        line = lines[line_idx].strip()
        if not line or line.startswith('!') or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f"Line {line_idx + 1} does not contain at least 3 k-point coordinates: {line}")
        
        # Take only first 3 columns (kx, ky, kz)
        kx, ky, kz = float(parts[0]), float(parts[1]), float(parts[2])
        kpoints.append([kx, ky, kz])
    
    kpoints = np.array(kpoints)
    
    if len(kpoints) != nk:
        raise ValueError(f"Expected {nk} k-points but found {len(kpoints)}")
    
    return kpoints, kpoints_type


def save_kpoints(kpoints, output_path, format='txt'):
    """
    Save k-points to file.
    
    Parameters
    ----------
    kpoints : np.ndarray
        Array of k-points
    output_path : str or Path
        Output file path
    format : str
        Output format: 'txt', 'numpy', or 'python'
    """
    output_path = Path(output_path)
    
    if format == 'txt':
        np.savetxt(output_path, kpoints, fmt='%16.8f')
        print(f"K-points saved to: {output_path}")
    
    elif format == 'numpy':
        np.save(output_path, kpoints)
        print(f"K-points saved to: {output_path}")
    
    elif format == 'python':
        with open(output_path, 'w') as f:
            f.write("kpoints = [\n")
            for kpt in kpoints:
                f.write(f"    [{kpt[0]:.8f}, {kpt[1]:.8f}, {kpt[2]:.8f}],\n")
            f.write("]\n")
        print(f"K-points saved to: {output_path}")
    
    else:
        raise ValueError(f"Unknown format: {format}")


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Extract k-points from Quantum ESPRESSO pw.x input files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract k-points to text file
  python extract_kpoints.py pwnscf.in -o kpoints.txt
  
  # Extract to numpy array
  python extract_kpoints.py pwnscf.in -o kpoints.npy -f numpy
  
  # Extract to Python list format
  python extract_kpoints.py pwnscf.in -o kpoints.py -f python
  
  # Print to stdout
  python extract_kpoints.py pwnscf.in
        """
    )
    
    parser.add_argument(
        'input_file',
        type=str,
        help="Path to pw.x input file (e.g., pwnscf.in)"
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help="Output file path (if not specified, prints to stdout)"
    )
    
    parser.add_argument(
        '-f', '--format',
        type=str,
        choices=['txt', 'numpy', 'python'],
        default='txt',
        help="Output format: txt (default), numpy (.npy), or python (list)"
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    try:
        # Extract k-points
        kpoints, kpoints_type = extract_kpoints_from_pwin(args.input_file)
        
        if args.verbose:
            print(f"K-points type: {kpoints_type}")
            print(f"Number of k-points: {len(kpoints)}")
        
        # Save or print
        if args.output:
            save_kpoints(kpoints, args.output, format=args.format)
        else:
            # Print to stdout
            print(f"# K-points extracted from: {args.input_file}")
            print(f"# Type: {kpoints_type}")
            print(f"# Number of k-points: {len(kpoints)}")
            for i, kpt in enumerate(kpoints):
                print(f"{kpt[0]:16.8f}  {kpt[1]:16.8f}  {kpt[2]:16.8f}")
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
