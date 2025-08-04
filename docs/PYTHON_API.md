# DNCS Python API Documentation

## Overview

The DNCS Python API provides a high-level interface for molecular conformational sampling and analysis. It includes functionality for polymer system creation, energy calculations, conformational sampling using Sobol sequences, and integration with OpenMM for molecular dynamics simulations.

## Installation

```bash
# Install Python dependencies
pip install -r requirements.txt

# Build and install the Python module
maturin develop --release
```

## Quick Start

```python
import dncs

# Create a polymer system
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# Calculate energy
energy = polymer.getEnergy()
print(f"Energy: {energy} kcal/mol")

# Generate conformational samples
sampler = dncs.SobolSampler(polymer, 100, "fold", 300.0, "output")

# Save structure
polymer.toPDB("structure.pdb")
```

## API Reference

### Functions

#### `getPDB(seq: str, filename: str) -> None`

Generates a PDB file from an amino acid sequence using the AmberFB15 force field.

**Parameters:**
- `seq` (str): Single-letter amino acid sequence (e.g., "YGGFM")
- `filename` (str): Output PDB filename

**Example:**
```python
import dncs

# Generate PDB for a pentapeptide
dncs.getPDB("YGGFM", "pentapeptide.pdb")

# Generate PDB for a longer sequence
dncs.getPDB("ALANYLGLYCYLPHENYLALANINE", "tripeptide.pdb")
```

**Supported Amino Acids:**
- A (Alanine), R (Arginine), N (Asparagine), D (Aspartic acid)
- C (Cysteine), Q (Glutamine), E (Glutamic acid), G (Glycine)
- H (Histidine), I (Isoleucine), L (Leucine), K (Lysine)
- M (Methionine), F (Phenylalanine), P (Proline), S (Serine)
- T (Threonine), W (Tryptophan), Y (Tyrosine), V (Valine)

#### `pdb_to_angle(filename: str) -> str`

Extracts dihedral angles from a PDB file and returns them as a comma-separated string.

**Parameters:**
- `filename` (str): Input PDB filename

**Returns:**
- `str`: Comma-separated dihedral angles in degrees

**Example:**
```python
# Extract dihedral angles from a PDB file
angles_csv = dncs.pdb_to_angle("structure.pdb")
print(f"Dihedral angles: {angles_csv}")

# Parse the angles for further processing
angles = [float(x.strip()) for x in angles_csv.split(",") if x.strip()]
print(f"Parsed angles: {angles}")
```

### Classes

#### `Polymer`

Represents a polymer system with force field parameters and provides methods for energy calculation and structure manipulation.

```python
class Polymer:
    def __init__(self, seq: str, forcefield: str) -> None
    def getEnergy(self) -> float
    def toPDB(self, filename: str) -> None
    def dihedral(self, foldername: str) -> None
```

##### Constructor

```python
Polymer(seq: str, forcefield: str) -> Polymer
```

Creates a new polymer system from an amino acid sequence.

**Parameters:**
- `seq` (str): Single-letter amino acid sequence
- `forcefield` (str): Force field name (must be one of the supported force fields)

**Supported Force Fields:**
- `"amber03.xml"` - AMBER03 force field
- `"amber10.xml"` - AMBER10 force field  
- `"amber96.xml"` - AMBER96 force field
- `"amber99sb.xml"` - AMBER99SB force field
- `"amberfb15.xml"` - AMBER-FB15 force field (recommended)

**Example:**
```python
# Create polymer with different force fields
polymer_fb15 = dncs.Polymer("YGGFM", "amberfb15.xml")
polymer_99sb = dncs.Polymer("ALANINE", "amber99sb.xml")

# Handle invalid force field
try:
    invalid_polymer = dncs.Polymer("YGG", "invalid.xml")
except SystemExit:
    print("Invalid force field specified")
```

##### Methods

###### `getEnergy() -> float`

Calculates the total potential energy of the polymer system.

**Returns:**
- `float`: Total energy in kcal/mol

**Example:**
```python
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")
energy = polymer.getEnergy()
print(f"System energy: {energy:.2f} kcal/mol")

# Compare energies of different conformations
angles1 = [0.0, 0.0, 0.0]  # Extended conformation
angles2 = [1.57, -1.57, 0.0]  # Partially folded

# Note: Direct angle setting not available in current API
# Use SobolSampler for conformation generation
```

###### `toPDB(filename: str) -> None`

Saves the current polymer structure to a PDB file.

**Parameters:**
- `filename` (str): Output PDB filename

**Example:**
```python
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# Save initial structure
polymer.toPDB("initial_structure.pdb")

# Save with timestamp
import datetime
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
polymer.toPDB(f"structure_{timestamp}.pdb")
```

###### `dihedral(foldername: str) -> None`

Logs dihedral angles to files in the specified folder.

**Parameters:**
- `foldername` (str): Output folder path

**Example:**
```python
import os

polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# Create output directory
os.makedirs("dihedral_output", exist_ok=True)

# Log dihedral angles
polymer.dihedral("dihedral_output")

# The folder will contain angle information files
```

#### `SobolSampler`

Performs conformational sampling using Sobol sequences for uniform exploration of conformational space.

```python
class SobolSampler:
    def __init__(self, system: Polymer, no_of_samples: int, method: str, temp: float, folder: str) -> None
```

##### Constructor

```python
SobolSampler(system: Polymer, no_of_samples: int, method: str, temp: float, folder: str) -> SobolSampler
```

Creates a new Sobol sampler for conformational sampling.

**Parameters:**
- `system` (Polymer): Polymer system to sample
- `no_of_samples` (int): Number of conformational samples to generate
- `method` (str): Sampling method ("fold", "search", or "explore")
- `temp` (float): Temperature in Kelvin for Boltzmann weighting
- `folder` (str): Output directory for sampled structures

**Sampling Methods:**
- `"fold"` - Optimized for protein folding simulations
- `"search"` - General conformational search
- `"explore"` - Extensive exploration of conformational space

**Example:**
```python
import os

# Create polymer system
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# Create output directory
os.makedirs("sampling_output", exist_ok=True)

# Create sampler with different methods
fold_sampler = dncs.SobolSampler(
    polymer, 1000, "fold", 300.0, "sampling_output/fold"
)

search_sampler = dncs.SobolSampler(
    polymer, 500, "search", 310.0, "sampling_output/search"
)

explore_sampler = dncs.SobolSampler(
    polymer, 2000, "explore", 298.0, "sampling_output/explore"
)
```

## Integration with OpenMM

The DNCS Python API can be integrated with OpenMM for molecular dynamics simulations through the `integrator.py` module.

### Configuration

Create a `dncs.toml` configuration file:

```toml
[simulation]
moleculename = "MyProtein"
folder = "Results"
sequence = "YGGFM"
interface = "openmm"
n_samples = 100
md_simulation = 10
temp = 300.0
forcefield = ["amber14-all.xml", "amber14/tip3p.xml"]
device = "CUDA"  # or "CPU", "HIP"
solvent = 10
steps = 50
gamma = 1.0
dt = 0.002
md_steps = 1
method = "fold"
```

### Usage with OpenMM

```python
import toml
import dncs
from integrator import DncsIntegrator

# Load configuration
config_data = toml.load("dncs.toml")
config = type('Config', (), config_data["simulation"])()

# Generate initial samples using DNCS
polymer = dncs.Polymer(config.sequence, "amberfb15.xml")
sampler = dncs.SobolSampler(
    polymer,
    config.n_samples,
    config.method,
    config.temp,
    f"{config.folder}/{config.moleculename}"
)

# Run OpenMM integration
integrator = DncsIntegrator(config)
integrator.run_integrator()
```

## Complete Examples

### Basic Conformational Analysis

```python
import dncs
import os
import numpy as np

def analyze_conformation(sequence, forcefield="amberfb15.xml"):
    """Analyze conformational properties of a peptide."""
    
    # Create polymer system
    polymer = dncs.Polymer(sequence, forcefield)
    
    # Calculate initial energy
    initial_energy = polymer.getEnergy()
    print(f"Initial energy: {initial_energy:.2f} kcal/mol")
    
    # Save initial structure
    polymer.toPDB("initial.pdb")
    
    # Extract initial dihedral angles
    initial_angles = dncs.pdb_to_angle("initial.pdb")
    print(f"Initial angles: {initial_angles}")
    
    # Create output directory
    output_dir = f"analysis_{sequence}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Perform conformational sampling
    sampler = dncs.SobolSampler(
        polymer, 500, "fold", 300.0, output_dir
    )
    
    print(f"Conformational sampling completed in {output_dir}/")
    
    return output_dir

# Analyze different peptides
sequences = ["YGG", "YGGFM", "ALANINE"]
for seq in sequences:
    print(f"\nAnalyzing {seq}:")
    analyze_conformation(seq)
```

### Temperature-Dependent Sampling

```python
import dncs
import os

def temperature_study(sequence, temperatures, samples_per_temp=200):
    """Study conformational behavior at different temperatures."""
    
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    results = {}
    
    for temp in temperatures:
        print(f"Sampling at {temp}K...")
        
        output_dir = f"temp_study_{sequence}_{temp}K"
        os.makedirs(output_dir, exist_ok=True)
        
        sampler = dncs.SobolSampler(
            polymer, samples_per_temp, "fold", temp, output_dir
        )
        
        results[temp] = output_dir
        print(f"Completed sampling at {temp}K")
    
    return results

# Study temperature effects
temperatures = [250.0, 300.0, 350.0, 400.0]
results = temperature_study("YGGFM", temperatures)

print("\nTemperature study results:")
for temp, folder in results.items():
    print(f"{temp}K: {folder}")
```

### Force Field Comparison

```python
import dncs
import os

def compare_forcefields(sequence, forcefields):
    """Compare different force fields for the same sequence."""
    
    results = {}
    
    for ff in forcefields:
        print(f"Testing {ff}...")
        
        try:
            polymer = dncs.Polymer(sequence, ff)
            energy = polymer.getEnergy()
            
            # Save structure
            ff_name = ff.replace(".xml", "")
            filename = f"{sequence}_{ff_name}.pdb"
            polymer.toPDB(filename)
            
            # Get dihedral angles
            angles = dncs.pdb_to_angle(filename)
            
            results[ff] = {
                'energy': energy,
                'angles': angles,
                'pdb_file': filename
            }
            
            print(f"{ff}: Energy = {energy:.2f} kcal/mol")
            
        except SystemExit:
            print(f"{ff}: Not supported or error occurred")
            results[ff] = None
    
    return results

# Compare force fields
forcefields = [
    "amber03.xml",
    "amber10.xml", 
    "amber96.xml",
    "amber99sb.xml",
    "amberfb15.xml"
]

comparison = compare_forcefields("YGG", forcefields)

print("\nForce field comparison:")
for ff, data in comparison.items():
    if data:
        print(f"{ff}: {data['energy']:.2f} kcal/mol")
    else:
        print(f"{ff}: Failed")
```

### Batch Processing

```python
import dncs
import os
from concurrent.futures import ThreadPoolExecutor
import time

def process_sequence(seq_data):
    """Process a single sequence."""
    sequence, output_base = seq_data
    
    try:
        # Create polymer
        polymer = dncs.Polymer(sequence, "amberfb15.xml")
        
        # Calculate energy
        energy = polymer.getEnergy()
        
        # Create output directory
        output_dir = f"{output_base}/{sequence}"
        os.makedirs(output_dir, exist_ok=True)
        
        # Save initial structure
        polymer.toPDB(f"{output_dir}/initial.pdb")
        
        # Perform sampling
        sampler = dncs.SobolSampler(
            polymer, 100, "fold", 300.0, output_dir
        )
        
        return {
            'sequence': sequence,
            'energy': energy,
            'output_dir': output_dir,
            'status': 'success'
        }
        
    except Exception as e:
        return {
            'sequence': sequence,
            'error': str(e),
            'status': 'failed'
        }

def batch_process_sequences(sequences, output_base="batch_output", max_workers=4):
    """Process multiple sequences in parallel."""
    
    os.makedirs(output_base, exist_ok=True)
    
    # Prepare sequence data
    seq_data = [(seq, output_base) for seq in sequences]
    
    results = []
    start_time = time.time()
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = executor.map(process_sequence, seq_data)
        results = list(futures)
    
    end_time = time.time()
    
    # Print summary
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']
    
    print(f"\nBatch processing completed in {end_time - start_time:.2f} seconds")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    
    if successful:
        print("\nSuccessful sequences:")
        for result in successful:
            print(f"  {result['sequence']}: {result['energy']:.2f} kcal/mol")
    
    if failed:
        print("\nFailed sequences:")
        for result in failed:
            print(f"  {result['sequence']}: {result['error']}")
    
    return results

# Process multiple sequences
sequences = [
    "YGG", "YGGFM", "ALANINE", "GLYCINE", 
    "PHENYLALANINE", "TRYPTOPHAN", "ARGININE"
]

batch_results = batch_process_sequences(sequences)
```

## Error Handling

### Common Errors and Solutions

1. **Invalid Force Field**
```python
try:
    polymer = dncs.Polymer("YGG", "invalid_forcefield.xml")
except SystemExit:
    print("Error: Unsupported force field")
    # Use a supported force field instead
    polymer = dncs.Polymer("YGG", "amberfb15.xml")
```

2. **Invalid Sequence**
```python
# DNCS expects single-letter amino acid codes
valid_sequence = "YGGFM"      # Correct
invalid_sequence = "YGG-FM"   # Incorrect (contains hyphen)

try:
    polymer = dncs.Polymer(invalid_sequence, "amberfb15.xml")
except Exception as e:
    print(f"Sequence error: {e}")
```

3. **File I/O Errors**
```python
import os

# Ensure output directory exists
output_dir = "my_output"
os.makedirs(output_dir, exist_ok=True)

# Check file permissions before writing
try:
    polymer = dncs.Polymer("YGG", "amberfb15.xml")
    polymer.toPDB(f"{output_dir}/structure.pdb")
except PermissionError:
    print("Error: No write permission for output directory")
```

## Performance Tips

1. **Use appropriate sample sizes**
```python
# For quick testing
quick_sampler = dncs.SobolSampler(polymer, 50, "fold", 300.0, "test")

# For production runs
production_sampler = dncs.SobolSampler(polymer, 1000, "fold", 300.0, "production")
```

2. **Choose optimal sampling methods**
```python
# For protein folding studies
fold_sampler = dncs.SobolSampler(polymer, 1000, "fold", 300.0, "fold_study")

# For general conformational search
search_sampler = dncs.SobolSampler(polymer, 500, "search", 300.0, "search_study")
```

3. **Parallel processing for multiple sequences**
```python
# Use ThreadPoolExecutor for I/O-bound operations
# Use ProcessPoolExecutor for CPU-bound operations
from concurrent.futures import ProcessPoolExecutor

def process_in_parallel(sequences):
    with ProcessPoolExecutor() as executor:
        results = executor.map(analyze_sequence, sequences)
    return list(results)
```

## Integration Examples

### With NumPy and SciPy

```python
import dncs
import numpy as np
from scipy import stats

def statistical_analysis(sequence, n_samples=500):
    """Perform statistical analysis of conformational sampling."""
    
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    
    # Create sampler
    sampler = dncs.SobolSampler(
        polymer, n_samples, "fold", 300.0, f"stats_{sequence}"
    )
    
    # Analyze initial structure
    polymer.toPDB("temp_structure.pdb")
    angles_str = dncs.pdb_to_angle("temp_structure.pdb")
    
    if angles_str.strip():
        angles = np.array([float(x) for x in angles_str.split(",") if x.strip()])
        
        print(f"Sequence: {sequence}")
        print(f"Number of dihedral angles: {len(angles)}")
        print(f"Mean angle: {np.mean(angles):.2f}°")
        print(f"Std deviation: {np.std(angles):.2f}°")
        print(f"Range: {np.min(angles):.2f}° to {np.max(angles):.2f}°")
        
        return {
            'sequence': sequence,
            'angles': angles,
            'mean': np.mean(angles),
            'std': np.std(angles),
            'range': (np.min(angles), np.max(angles))
        }
    
    return None

# Analyze multiple sequences
sequences = ["YGG", "YGGFM", "ALANINE"]
for seq in sequences:
    result = statistical_analysis(seq)
    if result:
        print(f"\nStatistics for {seq}:")
        print(f"  Mean: {result['mean']:.2f}°")
        print(f"  Std:  {result['std']:.2f}°")
```

### With Matplotlib for Visualization

```python
import dncs
import matplotlib.pyplot as plt
import numpy as np

def plot_energy_comparison(sequences, forcefield="amberfb15.xml"):
    """Plot energy comparison for different sequences."""
    
    energies = []
    labels = []
    
    for seq in sequences:
        try:
            polymer = dncs.Polymer(seq, forcefield)
            energy = polymer.getEnergy()
            energies.append(energy)
            labels.append(seq)
        except:
            print(f"Skipping {seq} due to error")
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, energies, color='skyblue', edgecolor='navy', alpha=0.7)
    
    # Add value labels on bars
    for bar, energy in zip(bars, energies):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f'{energy:.1f}', ha='center', va='bottom')
    
    plt.xlabel('Amino Acid Sequence')
    plt.ylabel('Energy (kcal/mol)')
    plt.title(f'Energy Comparison ({forcefield})')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    
    return dict(zip(labels, energies))

# Plot comparison
sequences = ["YGG", "YGGFM", "ALANINE", "GLYCINE"]
energy_data = plot_energy_comparison(sequences)
```

This comprehensive Python API documentation provides detailed information about all available functions and classes, along with practical examples and best practices for using the DNCS Python interface.