#!/usr/bin/env python3
"""
PyMOL Script Generator for DNCS Sample Visualization

This script reads the dncs.log file from a simulation folder and creates a PyMOL script
to load and visualize the top N samples based on their minimized energies.

Usage:
    python scripts/create_pymol_script.py <FOLDER_NAME> <TOP_N>

Example:
    python scripts/create_pymol_script.py 6RRO 5

This will create a pymol.pml file that loads the 5 best (lowest energy) structures
from the Result/6RRO/sample/ directory.
"""

import sys
import os
import re
from pathlib import Path
from typing import List, Tuple

def parse_log_file(log_path: str) -> List[Tuple[int, float]]:
    """
    Parse the dncs.log file to extract minimized energies for each model.

    Returns:
        List of tuples (model_number, minimized_energy)
    """
    minimized_energies = []

    try:
        with open(log_path, 'r') as f:
            content = f.read()

        # Look for lines like "MINIMIZED ENERGY FOR MODEL X = Y kJ/mol"
        pattern = r'MINIMIZED ENERGY FOR MODEL (\d+) = ([-\d\.e\+]+) kJ/mol'
        matches = re.findall(pattern, content)

        for match in matches:
            model_num = int(match[0])
            energy = float(match[1])
            minimized_energies.append((model_num, energy))

    except FileNotFoundError:
        print(f"Error: Log file not found at {log_path}")
        return []
    except Exception as e:
        print(f"Error parsing log file: {e}")
        return []

    return minimized_energies

def find_sample_files(sample_dir: str) -> List[str]:
    """
    Find all sample_XXXX.pdb files in the sample directory.

    Returns:
        List of sample file paths that actually exist
    """
    sample_files = []
    sample_path = Path(sample_dir)

    if not sample_path.exists():
        print(f"Warning: Sample directory {sample_dir} does not exist")
        return []

    # Look for files matching pattern sample_XXXX.pdb
    for pdb_file in sample_path.glob("sample_*.pdb"):
        sample_files.append(str(pdb_file))

    return sorted(sample_files)

def create_pymol_script(folder_name: str, top_n: int) -> bool:
    """
    Create a PyMOL script to visualize the top N samples.

    Args:
        folder_name: Name of the folder in Result/ directory
        top_n: Number of top samples to visualize

    Returns:
        True if successful, False otherwise
    """
    # Paths
    result_dir = Path("Result") / folder_name
    log_path = result_dir / "dncs.log"
    sample_dir = result_dir / "sample"

    # Check if directories exist
    if not result_dir.exists():
        print(f"Error: Result directory {result_dir} does not exist")
        return False

    # Parse log file for energies
    print(f"Parsing log file: {log_path}")
    energies = parse_log_file(str(log_path))

    if not energies:
        print("No energy data found in log file")
        return False

    # Sort by energy (lowest first)
    energies.sort(key=lambda x: x[1])

    print(f"Found {len(energies)} structures with energy data")
    print("Top 10 energies:")
    for i, (model, energy) in enumerate(energies[:10]):
        print(f"  {i+1}. Model {model}: {energy:.2f} kJ/mol")

    # Get top N structures
    top_structures = energies[:top_n]

    # Find available sample files
    sample_files = find_sample_files(str(sample_dir))

    # Check for sample files and provide alternatives if not found
    available_models = []
    missing_samples = []

    # First, try to find sample files based on model numbers
    for model_num, energy in top_structures:
        sample_file = sample_dir / f"sample_{model_num:04d}.pdb"
        if sample_file.exists():
            available_models.append((model_num, str(sample_file)))
        else:
            missing_samples.append(model_num)

    print(f"Found {len(available_models)} out of {len(top_structures)} requested sample files")

    if missing_samples:
        print(f"Missing sample files for models: {missing_samples}")

        # Check for alternative files in the result directory
        alternatives_found = []

        # Check for minimized files
        minimized_dir = result_dir / "Minimized"
        if minimized_dir.exists():
            for model_num in missing_samples:
                minimized_file = minimized_dir / f"Minimized_{model_num:04d}.pdb"
                if minimized_file.exists():
                    available_models.append((model_num, str(minimized_file)))
                    alternatives_found.append(f"Model {model_num}: Using minimized structure")

        # Check for equilibrated files
        langevin_dir = result_dir / "Langevin"
        if langevin_dir.exists() and len(alternatives_found) < len(missing_samples):
            remaining_missing = [m for m in missing_samples if not any(f"Model {m}" in alt for alt in alternatives_found)]
            for model_num in remaining_missing:
                equilibrated_file = langevin_dir / f"Equilibrated_{model_num:04d}.pdb"
                if equilibrated_file.exists():
                    available_models.append((model_num, str(equilibrated_file)))
                    alternatives_found.append(f"Model {model_num}: Using equilibrated structure")

        if alternatives_found:
            print("Using alternative files:")
            for alt in alternatives_found:
                print(f"  ✓ {alt}")

    if len(sample_files) > 0:
        print(f"Sample directory contains {len(sample_files)} total sample files")

    if not available_models:
        print("Error: No sample, minimized, or equilibrated files found for top structures")
        print("\nAvailable options:")
        print("1. Check if sampling was completed successfully")
        print("2. Try a different folder name")
        print("3. Check if files exist in subdirectories:")
        for subdir in ["sample", "Minimized", "Langevin"]:
            subdir_path = result_dir / subdir
            if subdir_path.exists():
                pdb_files = list(subdir_path.glob("*.pdb"))
                print(f"   - {subdir}/: {len(pdb_files)} PDB files")
        return False

    # Create PyMOL script - simplified to just load models
    pymol_script = f"""# PyMOL Script for DNCS Sample Visualization
# Generated for folder: {folder_name}
# Top {top_n} structures by minimized energy

# Clear any existing objects
delete all

# Load structures
"""

    for i, (model_num, file_path) in enumerate(available_models):
        # Convert absolute path to relative path from root directory
        rel_path = os.path.relpath(file_path, ".")
        pymol_script += f"load {rel_path}, model_{model_num}\n"

    # Write the script
    script_path = "pymol.pml"
    try:
        with open(script_path, 'w') as f:
            f.write(pymol_script)

        print(f"\nPyMOL script created: {script_path}")
        print(f"Loaded {len(available_models)} structures from top {top_n} requested")
        print("\nTo use:")
        print("  pymol pymol.pml")
        print("  or")
        print("  pymol -c -d 'run pymol.pml'")

        return True

    except Exception as e:
        print(f"Error writing PyMOL script: {e}")
        return False

def main():
    if len(sys.argv) != 3:
        print("Usage: python scripts/create_pymol_script.py <FOLDER_NAME> <TOP_N>")
        print("\nExample:")
        print("  python scripts/create_pymol_script.py 6RRO 5")
        print("\nThis will create a pymol.pml script to load the 5 best structures")
        print("from Result/6RRO/sample/ based on minimized energies in dncs.log")
        sys.exit(1)

    folder_name = sys.argv[1]
    try:
        top_n = int(sys.argv[2])
        if top_n <= 0:
            raise ValueError("TOP_N must be positive")
    except ValueError as e:
        print(f"Error: Invalid TOP_N value '{sys.argv[2]}': {e}")
        sys.exit(1)

    print(f"Creating PyMOL visualization for top {top_n} samples from {folder_name}")

    success = create_pymol_script(folder_name, top_n)

    if success:
        sys.exit(0)
    else:
        print("\nFailed to create PyMOL script")
        sys.exit(1)

if __name__ == "__main__":
    main()
