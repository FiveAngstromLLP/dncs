#!/usr/bin/env python3
"""
DNCS Result Explorer

This script helps explore the available simulation results, showing what files
exist in each result folder and providing energy summaries.

Usage:
    python scripts/list_results.py                 # List all result folders
    python scripts/list_results.py <FOLDER_NAME>   # Detailed view of specific folder

Examples:
    python scripts/list_results.py                 # Show all available results
    python scripts/list_results.py 6RRO            # Detailed view of 6RRO results
"""

import sys
import os
import re
from pathlib import Path
from typing import List, Tuple, Dict

def parse_energies_from_log(log_path: str) -> Dict[str, List[Tuple[int, float]]]:
    """
    Parse energies from dncs.log file.

    Returns:
        Dictionary with keys: 'initial', 'minimized', 'equilibrated'
        Each value is a list of (model_number, energy) tuples
    """
    energies = {
        'initial': [],
        'minimized': [],
        'equilibrated': []
    }

    if not os.path.exists(log_path):
        return energies

    try:
        with open(log_path, 'r') as f:
            content = f.read()

        # Parse different energy types
        patterns = {
            'initial': r'ENERGY FOR MODEL (\d+) = ([-\d\.e\+]+) kJ/mol',
            'minimized': r'MINIMIZED ENERGY FOR MODEL (\d+) = ([-\d\.e\+]+) kJ/mol',
            'equilibrated': r'EQUILIBRATED ENERGY AFTER \d+ STEPS FOR MODEL (\d+) = ([-\d\.e\+]+) kJ/mol'
        }

        for energy_type, pattern in patterns.items():
            matches = re.findall(pattern, content)
            for match in matches:
                model_num = int(match[0])
                energy = float(match[1])
                energies[energy_type].append((model_num, energy))

    except Exception as e:
        print(f"Warning: Error parsing log file {log_path}: {e}")

    return energies

def count_files_in_directory(dir_path: Path, pattern: str = "*.pdb") -> int:
    """Count files matching pattern in directory."""
    if not dir_path.exists():
        return 0
    try:
        return len(list(dir_path.glob(pattern)))
    except:
        return 0

def get_folder_summary(folder_path: Path) -> Dict:
    """Get summary information about a result folder."""
    summary = {
        'name': folder_path.name,
        'exists': folder_path.exists(),
        'log_exists': False,
        'energies': {},
        'file_counts': {},
        'total_size_mb': 0
    }

    if not folder_path.exists():
        return summary

    # Check for log file
    log_path = folder_path / "dncs.log"
    summary['log_exists'] = log_path.exists()

    if summary['log_exists']:
        summary['energies'] = parse_energies_from_log(str(log_path))

    # Count files in subdirectories
    subdirs = ['sample', 'Minimized', 'Langevin', 'MDSimulation']
    for subdir in subdirs:
        subdir_path = folder_path / subdir
        summary['file_counts'][subdir] = count_files_in_directory(subdir_path)

    # Calculate total size
    try:
        total_size = sum(f.stat().st_size for f in folder_path.rglob('*') if f.is_file())
        summary['total_size_mb'] = total_size / (1024 * 1024)
    except:
        summary['total_size_mb'] = 0

    return summary

def list_all_results():
    """List all available result folders."""
    result_dir = Path("Result")

    if not result_dir.exists():
        print("Error: Result directory does not exist")
        print("   Run a simulation first to generate results")
        return

    folders = [f for f in result_dir.iterdir() if f.is_dir()]

    if not folders:
        print("Result directory is empty")
        print("   Run a simulation first to generate results")
        return

    print("Available DNCS Result Folders:")
    print("=" * 60)

    for folder in sorted(folders):
        summary = get_folder_summary(folder)

        if summary['log_exists']:
            # Show energy summary
            initial_count = len(summary['energies']['initial'])
            minimized_count = len(summary['energies']['minimized'])
            equilibrated_count = len(summary['energies']['equilibrated'])

            print(f"   Energies: {initial_count} initial, {minimized_count} minimized, {equilibrated_count} equilibrated")

            # Show best energy
            if minimized_count > 0:
                best_energy = min(summary['energies']['minimized'], key=lambda x: x[1])
                print(f"   Best energy: Model {best_energy[0]} = {best_energy[1]:.2f} kJ/mol")

        # Show file counts
        file_info = []
        for subdir, count in summary['file_counts'].items():
            if count > 0:
                file_info.append(f"{subdir}: {count}")

        if file_info:
            print(f"   Files: {', '.join(file_info)}")

        if summary['total_size_mb'] > 0:
            print(f"   Size: {summary['total_size_mb']:.1f} MB")

    print("\n" + "=" * 60)
    print("Usage:")
    print("   just list <FOLDER_NAME>     - Detailed view of specific folder")
    print("   just view <FOLDER_NAME> N   - Create PyMOL script for top N structures")

def show_detailed_folder_info(folder_name: str):
    """Show detailed information about a specific folder."""
    folder_path = Path("Result") / folder_name

    if not folder_path.exists():
        print(f"Error: Folder '{folder_name}' does not exist in Result/")
        list_all_results()
        return

    summary = get_folder_summary(folder_path)

    print(f"Detailed View: {folder_name}")
    print("=" * 60)

    # Basic info
    print(f"Path: Result/{folder_name}")
    print(f"Size: {summary['total_size_mb']:.1f} MB")
    print(f"Log file: {'Found' if summary['log_exists'] else 'Missing'}")

    if not summary['log_exists']:
        print("   Warning: No dncs.log found - simulation may not have completed")
        return

    # Energy analysis
    print("\nEnergy Analysis:")

    energies = summary['energies']
    for energy_type in ['initial', 'minimized', 'equilibrated']:
        data = energies[energy_type]
        if data:
            sorted_data = sorted(data, key=lambda x: x[1])
            best = sorted_data[0]
            worst = sorted_data[-1]

            print(f"   {energy_type.title()}:")
            print(f"     Count: {len(data)}")
            print(f"     Best:  Model {best[0]} = {best[1]:.2f} kJ/mol")
            print(f"     Worst: Model {worst[0]} = {worst[1]:.2f} kJ/mol")

            # Show top 5
            if len(sorted_data) > 1:
                print("     Top 5: ", end="")
                top5 = sorted_data[:5]
                print(", ".join([f"M{m}({e:.0f})" for m, e in top5]))

    # File structure
    print("\nFile Structure:")
    subdirs = {
        'sample': 'Initial sample structures',
        'Minimized': 'Energy minimized structures',
        'Langevin': 'Equilibrated structures',
        'MDSimulation': 'Production MD results'
    }

    for subdir, description in subdirs.items():
        count = summary['file_counts'][subdir]
        status = "[OK]" if count > 0 else "[EMPTY]"
        print(f"   {status} {subdir}/: {count} files - {description}")

    # Show available sample files for visualization
    sample_count = summary['file_counts']['sample']
    minimized_count = summary['file_counts']['Minimized']

    print("\nVisualization Options:")

    if sample_count > 0:
        print(f"   Sample files: {sample_count} available")
        print(f"      Command: just view {folder_name} <N>")
    elif minimized_count > 0:
        print(f"   Minimized files: {minimized_count} available")
        print(f"      Command: just view {folder_name} <N>  (will use minimized structures)")
    else:
        print("   No visualization files found")

    # Recommendations
    print("\nRecommendations:")

    if energies['minimized']:
        best_models = sorted(energies['minimized'], key=lambda x: x[1])[:5]
        model_numbers = [str(m) for m, e in best_models]
        print(f"   Best models to visualize: {', '.join(model_numbers)}")
        print(f"   Command: just view {folder_name} 5")

    # Check for issues
    issues = []
    if len(energies['initial']) != len(energies['minimized']):
        issues.append("Some structures failed minimization")
    if len(energies['minimized']) != len(energies['equilibrated']):
        issues.append("Some structures failed equilibration")

    if issues:
        print("   Issues detected:")
        for issue in issues:
            print(f"      - {issue}")

def main():
    if len(sys.argv) == 1:
        # No arguments - list all folders
        list_all_results()
    elif len(sys.argv) == 2:
        # One argument - show detailed folder info
        folder_name = sys.argv[1]
        show_detailed_folder_info(folder_name)
    else:
        print("Usage:")
        print("  python scripts/list_results.py                 # List all result folders")
        print("  python scripts/list_results.py <FOLDER_NAME>   # Detailed view of folder")
        print("\nExamples:")
        print("  python scripts/list_results.py")
        print("  python scripts/list_results.py 6RRO")
        sys.exit(1)

if __name__ == "__main__":
    main()
