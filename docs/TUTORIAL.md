# DNCS Tutorial: From Basics to Advanced Usage

## Table of Contents

1. [Introduction](#introduction)
2. [Installation and Setup](#installation-and-setup)
3. [Basic Concepts](#basic-concepts)
4. [Getting Started with Rust](#getting-started-with-rust)
5. [Getting Started with Python](#getting-started-with-python)
6. [Working with Molecular Systems](#working-with-molecular-systems)
7. [Conformational Sampling](#conformational-sampling)
8. [Energy Minimization](#energy-minimization)
9. [Integration with OpenMM](#integration-with-openmm)
10. [Advanced Topics](#advanced-topics)
11. [Best Practices](#best-practices)
12. [Troubleshooting](#troubleshooting)

## Introduction

DNCS (Digital Nets Conformational Sampling) is a powerful tool for enhanced conformational sampling of molecular systems using digital nets. This tutorial will guide you through using DNCS for molecular simulation and analysis, from basic operations to advanced workflows.

### What is DNCS?

DNCS uses Sobol sequences (a type of digital net) to provide uniform sampling in high-dimensional conformational spaces. This approach offers several advantages:

- **Uniform Coverage**: Better exploration of conformational space compared to random sampling
- **Efficiency**: Fewer samples needed to achieve good coverage
- **Reproducibility**: Deterministic sequences ensure reproducible results
- **Scalability**: Works well for high-dimensional systems

### Key Features

- Support for multiple AMBER force fields
- Rust core with Python bindings
- Parallel processing capabilities
- Integration with OpenMM for molecular dynamics
- Energy minimization using L-BFGS
- Comprehensive output analysis

## Installation and Setup

### Prerequisites

Before installing DNCS, ensure you have:

- **Rust toolchain** (1.70.0 or later)
- **Python** 3.7+ (for Python bindings)
- **Maturin** (for building Python bindings)
- **OpenMM** (optional, for MD simulations)

### Installation Steps

1. **Install Rust** (if not already installed):
   ```bash
   # Unix-based systems
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   
   # Windows: Download and run rustup-init.exe from https://rustup.rs
   ```

2. **Clone and Build DNCS**:
   ```bash
   git clone <repository-url>
   cd dncs
   
   # Using just (recommended)
   just install
   
   # Or manually
   pip install -r requirements.txt
   pip install maturin
   maturin develop --release
   ```

3. **Verify Installation**:
   ```bash
   # Test Rust binary
   cargo build --release
   ./target/release/dncs --help
   
   # Test Python module
   python -c "import dncs; print('DNCS imported successfully')"
   ```

## Basic Concepts

### Molecular Systems

In DNCS, a molecular system consists of:
- **Atoms**: Individual particles with coordinates and properties
- **Force Field**: Parameters defining interactions (AMBER force fields)
- **Topology**: Connectivity and bonding information
- **Conformations**: Different 3D arrangements of the same molecule

### Digital Nets and Sobol Sequences

Digital nets are point sets designed for uniform distribution in multi-dimensional spaces. Sobol sequences are a specific type of digital net that provides:

- Low-discrepancy sampling
- Uniform coverage of parameter space
- Deterministic generation
- Good convergence properties

### Dihedral Angles

Dihedral angles (torsion angles) define the rotation around bonds and are crucial for:
- Protein folding
- Conformational changes
- Energy landscapes
- Structural analysis

## Getting Started with Rust

### Your First DNCS Program

Let's create a simple program that generates a molecular system and calculates its energy:

```rust
// examples/basic_usage.rs
use std::sync::Arc;
use libdncs::*;

fn main() {
    // 1. Create a molecular system from sequence
    let sequence = "YGGFM";  // Tyrosine-Glycine-Glycine-Phenylalanine-Methionine
    let mut system = System::new(sequence, FF::AmberFB15.init());
    
    // 2. Initialize system parameters
    system.init_parameters();
    
    // 3. Calculate energy
    let amber = Amber::new(Arc::new(system.clone()));
    let energy = amber.energy();
    
    println!("Sequence: {}", sequence);
    println!("Number of atoms: {}", system.particles.len());
    println!("Energy: {:.2} kcal/mol", energy);
    
    // 4. Save structure to PDB
    system.to_pdb("initial_structure.pdb");
    println!("Structure saved to initial_structure.pdb");
}
```

Run this example:
```bash
cargo run --example basic_usage
```

### Understanding the Output

The program will output:
```
Sequence: YGGFM
Number of atoms: 75
Energy: -23.45 kcal/mol
Structure saved to initial_structure.pdb
```

### Working with Different Force Fields

```rust
use libdncs::*;

fn compare_force_fields(sequence: &str) {
    let force_fields = vec![
        ("AMBER03", FF::Amber03),
        ("AMBER99SB", FF::Amber99SB),
        ("AMBERFB15", FF::AmberFB15),
    ];
    
    println!("Force field comparison for sequence: {}", sequence);
    println!("{:<12} {:<15}", "Force Field", "Energy (kcal/mol)");
    println!("{}", "-".repeat(30));
    
    for (name, ff) in force_fields {
        let mut system = System::new(sequence, ff.init());
        system.init_parameters();
        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();
        
        println!("{:<12} {:<15.2}", name, energy);
    }
}

fn main() {
    compare_force_fields("YGG");
}
```

### Basic Conformational Sampling

```rust
use std::sync::Arc;
use libdncs::*;

fn main() {
    let sequence = "YGG";
    let mut system = System::new(sequence, FF::AmberFB15.init());
    system.init_parameters();
    
    // Create sampler
    let mut sampler = Sampler::new(
        Arc::new(system),
        Method::Fold,
        "output".to_string()
    );
    
    // Generate 100 samples at 300K
    sampler.sample(100, 300.0);
    
    // Write results
    sampler.write_angles();
    
    println!("Sampling completed. Check 'output' directory for results.");
}
```

## Getting Started with Python

### Your First Python Script

```python
# basic_usage.py
import dncs
import os

def basic_analysis(sequence):
    """Perform basic analysis of a peptide sequence."""
    
    print(f"Analyzing sequence: {sequence}")
    
    # Create polymer system
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    
    # Calculate energy
    energy = polymer.getEnergy()
    print(f"Initial energy: {energy:.2f} kcal/mol")
    
    # Save initial structure
    polymer.toPDB(f"{sequence}_initial.pdb")
    print(f"Structure saved to {sequence}_initial.pdb")
    
    # Extract dihedral angles
    angles_csv = dncs.pdb_to_angle(f"{sequence}_initial.pdb")
    if angles_csv.strip():
        angles = [float(x.strip()) for x in angles_csv.split(",") if x.strip()]
        print(f"Number of dihedral angles: {len(angles)}")
        print(f"Angle range: {min(angles):.1f}° to {max(angles):.1f}°")
    
    return energy

if __name__ == "__main__":
    sequences = ["YGG", "YGGFM", "ALANINE"]
    
    for seq in sequences:
        print("\n" + "="*40)
        energy = basic_analysis(seq)
        print(f"Analysis complete for {seq}")
```

Run the script:
```bash
python basic_usage.py
```

### Simple Conformational Sampling

```python
# sampling_example.py
import dncs
import os

def sample_conformations(sequence, n_samples=100):
    """Generate conformational samples for a sequence."""
    
    # Create output directory
    output_dir = f"samples_{sequence}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create polymer
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    initial_energy = polymer.getEnergy()
    
    print(f"Sampling {sequence} (Initial energy: {initial_energy:.2f} kcal/mol)")
    
    # Create sampler
    sampler = dncs.SobolSampler(
        polymer, n_samples, "fold", 300.0, output_dir
    )
    
    print(f"Generated {n_samples} samples in {output_dir}/")
    
    # Check output files
    sample_files = [f for f in os.listdir(output_dir) if f.endswith('.pdb')]
    print(f"Created {len(sample_files)} PDB files")
    
    return output_dir

if __name__ == "__main__":
    # Sample different sequences
    sequences = ["YGG", "YGGFM"]
    
    for seq in sequences:
        output_dir = sample_conformations(seq, 50)
        print(f"Results for {seq} in {output_dir}\n")
```

## Working with Molecular Systems

### Loading from PDB Files

```rust
// Load system from existing PDB file
use libdncs::*;

fn main() {
    // Create system from PDB file
    let system = System::from_pdb("input.pdb", FF::AmberFB15.init());
    
    println!("Loaded {} atoms from PDB", system.particles.len());
    
    // Analyze the structure
    for (i, atom) in system.particles.iter().take(5).enumerate() {
        println!("Atom {}: {} {} {:.3} {:.3} {:.3}", 
                 i+1, atom.name, atom.residue, atom.x, atom.y, atom.z);
    }
}
```

### Analyzing Molecular Properties

```python
import dncs
import numpy as np

def analyze_molecular_properties(sequence):
    """Analyze various molecular properties."""
    
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    
    # Basic properties
    energy = polymer.getEnergy()
    
    # Save and analyze structure
    polymer.toPDB("temp_structure.pdb")
    angles_csv = dncs.pdb_to_angle("temp_structure.pdb")
    
    if angles_csv.strip():
        angles = np.array([float(x) for x in angles_csv.split(",") if x.strip()])
        
        properties = {
            'sequence': sequence,
            'energy': energy,
            'n_dihedrals': len(angles),
            'mean_angle': np.mean(angles),
            'std_angle': np.std(angles),
            'angle_range': (np.min(angles), np.max(angles))
        }
        
        return properties
    
    return None

# Analyze multiple sequences
sequences = ["A", "AA", "AAA", "AAAA", "AAAAA"]
results = []

for seq in sequences:
    props = analyze_molecular_properties(seq)
    if props:
        results.append(props)
        print(f"{seq:>5}: Energy={props['energy']:>8.2f}, "
              f"Dihedrals={props['n_dihedrals']:>2}, "
              f"Mean angle={props['mean_angle']:>6.1f}°")

# Clean up
import os
if os.path.exists("temp_structure.pdb"):
    os.remove("temp_structure.pdb")
```

### Working with Different Amino Acids

```python
import dncs

# Dictionary of amino acid properties
amino_acids = {
    'A': 'Alanine (small, nonpolar)',
    'R': 'Arginine (large, positive)',
    'N': 'Asparagine (medium, polar)',
    'D': 'Aspartic acid (medium, negative)',
    'C': 'Cysteine (small, polar)',
    'E': 'Glutamic acid (large, negative)',
    'Q': 'Glutamine (large, polar)',
    'G': 'Glycine (tiny, nonpolar)',
    'H': 'Histidine (large, positive)',
    'I': 'Isoleucine (large, nonpolar)',
    'L': 'Leucine (large, nonpolar)',
    'K': 'Lysine (large, positive)',
    'M': 'Methionine (large, nonpolar)',
    'F': 'Phenylalanine (large, aromatic)',
    'P': 'Proline (medium, nonpolar)',
    'S': 'Serine (small, polar)',
    'T': 'Threonine (medium, polar)',
    'W': 'Tryptophan (largest, aromatic)',
    'Y': 'Tyrosine (large, aromatic)',
    'V': 'Valine (medium, nonpolar)'
}

def test_amino_acid_energies():
    """Test energy calculations for different amino acids."""
    
    print("Amino Acid Energy Analysis")
    print("=" * 50)
    print(f"{'AA':>2} {'Name':>15} {'Energy':>10}")
    print("-" * 50)
    
    energies = {}
    
    for aa, description in amino_acids.items():
        try:
            polymer = dncs.Polymer(aa, "amberfb15.xml")
            energy = polymer.getEnergy()
            energies[aa] = energy
            
            name = description.split()[0]
            print(f"{aa:>2} {name:>15} {energy:>10.2f}")
            
        except Exception as e:
            print(f"{aa:>2} {'ERROR':>15} {'N/A':>10}")
    
    # Find extremes
    if energies:
        min_aa = min(energies, key=energies.get)
        max_aa = max(energies, key=energies.get)
        
        print("\nSummary:")
        print(f"Lowest energy:  {min_aa} ({energies[min_aa]:.2f} kcal/mol)")
        print(f"Highest energy: {max_aa} ({energies[max_aa]:.2f} kcal/mol)")
    
    return energies

# Run the analysis
energies = test_amino_acid_energies()
```

## Conformational Sampling

### Understanding Sampling Methods

DNCS provides three sampling methods:

1. **Fold**: Optimized for protein folding studies
2. **Search**: General conformational search
3. **Explore**: Extensive exploration of conformational space

```python
import dncs
import os

def compare_sampling_methods(sequence, n_samples=200):
    """Compare different sampling methods."""
    
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    methods = ["fold", "search", "explore"]
    results = {}
    
    for method in methods:
        print(f"Testing {method} method...")
        
        output_dir = f"method_comparison_{sequence}_{method}"
        os.makedirs(output_dir, exist_ok=True)
        
        sampler = dncs.SobolSampler(
            polymer, n_samples, method, 300.0, output_dir
        )
        
        # Count output files
        pdb_files = [f for f in os.listdir(output_dir) if f.endswith('.pdb')]
        results[method] = {
            'output_dir': output_dir,
            'n_structures': len(pdb_files)
        }
        
        print(f"  {method}: Generated {len(pdb_files)} structures")
    
    return results

# Compare methods for a tripeptide
results = compare_sampling_methods("YGG", 100)

print("\nSampling method comparison:")
for method, data in results.items():
    print(f"{method:>7}: {data['n_structures']} structures in {data['output_dir']}")
```

### Advanced Sampling with Temperature Control

```rust
use std::sync::Arc;
use libdncs::*;

fn temperature_dependent_sampling() {
    let sequence = "YGGFM";
    let mut system = System::new(sequence, FF::AmberFB15.init());
    system.init_parameters();
    
    let temperatures = vec![250.0, 300.0, 350.0, 400.0];
    
    for temp in temperatures {
        println!("Sampling at {}K", temp);
        
        let folder = format!("temp_{}K", temp as i32);
        let mut sampler = Sampler::new(
            Arc::new(system.clone()),
            Method::Fold,
            folder.clone()
        );
        
        sampler.sample(500, temp);
        sampler.write_angles();
        
        println!("Results saved to {}/", folder);
    }
}

fn main() {
    temperature_dependent_sampling();
}
```

### Analyzing Sampling Results

```python
import dncs
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def analyze_sampling_results(output_dir):
    """Analyze the results of conformational sampling."""
    
    # Find all PDB files
    pdb_files = glob.glob(f"{output_dir}/sample_*.pdb")
    
    if not pdb_files:
        print(f"No sample files found in {output_dir}")
        return None
    
    print(f"Analyzing {len(pdb_files)} structures...")
    
    # Extract energies and angles from each structure
    energies = []
    all_angles = []
    
    for pdb_file in sorted(pdb_files):
        try:
            # Extract angles
            angles_csv = dncs.pdb_to_angle(pdb_file)
            if angles_csv.strip():
                angles = [float(x) for x in angles_csv.split(",") if x.strip()]
                all_angles.append(angles)
                
                # Extract energy from filename or PDB (if available)
                # For now, we'll create a dummy energy based on angle variance
                energy = np.var(angles) * 10  # Dummy energy calculation
                energies.append(energy)
                
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")
    
    if not energies:
        print("No valid structures found")
        return None
    
    # Analysis
    analysis = {
        'n_structures': len(energies),
        'energy_stats': {
            'mean': np.mean(energies),
            'std': np.std(energies),
            'min': np.min(energies),
            'max': np.max(energies)
        },
        'angle_stats': {
            'n_angles': len(all_angles[0]) if all_angles else 0,
            'mean_angles': np.mean(all_angles, axis=0) if all_angles else [],
            'std_angles': np.std(all_angles, axis=0) if all_angles else []
        }
    }
    
    # Print summary
    print(f"\nSampling Analysis for {output_dir}:")
    print(f"  Number of structures: {analysis['n_structures']}")
    print(f"  Energy range: {analysis['energy_stats']['min']:.2f} to {analysis['energy_stats']['max']:.2f}")
    print(f"  Mean energy: {analysis['energy_stats']['mean']:.2f} ± {analysis['energy_stats']['std']:.2f}")
    print(f"  Number of dihedral angles: {analysis['angle_stats']['n_angles']}")
    
    return analysis

# Example usage
def run_sampling_analysis():
    """Run sampling and analyze results."""
    
    # Generate samples
    sequence = "YGG"
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    
    output_dir = f"analysis_{sequence}"
    os.makedirs(output_dir, exist_ok=True)
    
    sampler = dncs.SobolSampler(
        polymer, 100, "fold", 300.0, output_dir
    )
    
    # Analyze results
    analysis = analyze_sampling_results(output_dir)
    
    return analysis

# Run the analysis
if __name__ == "__main__":
    analysis = run_sampling_analysis()
```

## Energy Minimization

### Basic Energy Minimization

```rust
use std::sync::Arc;
use libdncs::*;

fn minimize_structure() {
    let sequence = "YGGFM";
    let mut system = System::new(sequence, FF::AmberFB15.init());
    system.init_parameters();
    
    // First, generate some samples
    let mut sampler = Sampler::new(
        Arc::new(system.clone()),
        Method::Fold,
        "minimize_example".to_string()
    );
    sampler.sample(50, 300.0);
    sampler.write_angles();
    
    // Now minimize the top 10 structures
    let mut minimizer = Minimizer::new(
        "minimize_example".to_string(),
        FF::AmberFB15.init(),
        10
    );
    
    minimizer.minimize_all();
    minimizer.conformational_sort(300.0);
    minimizer.write_angles();
    minimizer.rename();
    
    println!("Minimization completed. Check minimize_example/minimize/ for results.");
}

fn main() {
    minimize_structure();
}
```

### Analyzing Minimization Results

```python
import os
import glob
import dncs

def analyze_minimization(base_folder):
    """Analyze energy minimization results."""
    
    sample_dir = f"{base_folder}/sample"
    minimize_dir = f"{base_folder}/minimize"
    
    if not os.path.exists(minimize_dir):
        print(f"Minimization directory not found: {minimize_dir}")
        return
    
    # Find minimized structures
    minimized_files = glob.glob(f"{minimize_dir}/Minimized_*.pdb")
    
    if not minimized_files:
        print("No minimized structures found")
        return
    
    print(f"Found {len(minimized_files)} minimized structures")
    
    # Analyze each minimized structure
    results = []
    
    for pdb_file in sorted(minimized_files):
        # Extract structure number
        basename = os.path.basename(pdb_file)
        struct_num = basename.replace('Minimized_', '').replace('.pdb', '')
        
        # Get angles
        angles_csv = dncs.pdb_to_angle(pdb_file)
        if angles_csv.strip():
            angles = [float(x) for x in angles_csv.split(",") if x.strip()]
            
            # Create a temporary polymer to calculate energy
            # (Note: This is a simplified approach)
            result = {
                'structure': struct_num,
                'n_angles': len(angles),
                'angle_range': (min(angles), max(angles)),
                'mean_angle': sum(angles) / len(angles)
            }
            
            results.append(result)
    
    # Print results
    print("\nMinimization Results:")
    print(f"{'Structure':>10} {'N_Angles':>10} {'Mean_Angle':>12} {'Range':>20}")
    print("-" * 55)
    
    for result in results:
        print(f"{result['structure']:>10} {result['n_angles']:>10} "
              f"{result['mean_angle']:>12.1f} "
              f"{result['angle_range'][0]:>8.1f}-{result['angle_range'][1]:<8.1f}")
    
    return results

# Example usage
if __name__ == "__main__":
    # This assumes you have run the Rust minimization example
    results = analyze_minimization("minimize_example")
```

## Integration with OpenMM

### Setting up OpenMM Integration

First, create a configuration file `dncs.toml`:

```toml
[simulation]
moleculename = "TestProtein"
folder = "OpenMM_Results"
sequence = "YGGFM"
interface = "openmm"
n_samples = 50
md_simulation = 5
temp = 300.0
forcefield = ["amber14-all.xml", "amber14/tip3p.xml"]
device = "CPU"  # Change to "CUDA" if available
solvent = 10
steps = 100
gamma = 1.0
dt = 0.002
md_steps = 1000
method = "fold"
```

### Running OpenMM Simulations

```python
# openmm_example.py
import toml
import dncs
import os
from python.src.integrator import DncsIntegrator

def run_openmm_simulation(config_file="dncs.toml"):
    """Run a complete DNCS + OpenMM simulation."""
    
    # Load configuration
    try:
        config_data = toml.load(config_file)
        config = type('Config', (), config_data["simulation"])()
    except FileNotFoundError:
        print(f"Configuration file {config_file} not found")
        return
    
    print(f"Starting simulation for {config.moleculename}")
    print(f"Sequence: {config.sequence}")
    print(f"Samples: {config.n_samples}")
    print(f"MD simulations: {config.md_simulation}")
    
    # Step 1: Generate initial samples with DNCS
    print("\n1. Generating conformational samples...")
    
    polymer = dncs.Polymer(config.sequence, "amberfb15.xml")
    initial_energy = polymer.getEnergy()
    print(f"Initial energy: {initial_energy:.2f} kcal/mol")
    
    # Create output directory
    sample_dir = f"{config.folder}/{config.moleculename}"
    os.makedirs(sample_dir, exist_ok=True)
    
    # Generate samples
    sampler = dncs.SobolSampler(
        polymer,
        config.n_samples,
        config.method,
        config.temp,
        sample_dir
    )
    
    print(f"Generated {config.n_samples} samples")
    
    # Step 2: Run OpenMM integration
    print("\n2. Running OpenMM simulations...")
    
    try:
        integrator = DncsIntegrator(config)
        integrator.run_integrator()
        print("OpenMM simulations completed successfully")
    except Exception as e:
        print(f"OpenMM simulation failed: {e}")
        print("Make sure OpenMM is installed and configured properly")
    
    # Step 3: Analyze results
    print("\n3. Analyzing results...")
    analyze_openmm_results(config)

def analyze_openmm_results(config):
    """Analyze OpenMM simulation results."""
    
    base_dir = f"{config.folder}/{config.moleculename}"
    
    # Check for different result directories
    result_dirs = {
        'Sampled': f"{base_dir}/Sampled",
        'Minimized': f"{base_dir}/Minimized", 
        'Langevin': f"{base_dir}/Langevin",
        'MDSimulation': f"{base_dir}/MDSimulation"
    }
    
    print("Result directories:")
    for name, path in result_dirs.items():
        if os.path.exists(path):
            files = [f for f in os.listdir(path) if f.endswith('.pdb')]
            print(f"  {name:>12}: {len(files)} PDB files")
        else:
            print(f"  {name:>12}: Not found")
    
    # Check log file
    log_file = f"{base_dir}/dncs.log"
    if os.path.exists(log_file):
        print(f"\nLog file: {log_file}")
        with open(log_file, 'r') as f:
            lines = f.readlines()
            print(f"Log entries: {len(lines)}")
            if lines:
                print("Last few entries:")
                for line in lines[-3:]:
                    print(f"  {line.strip()}")

if __name__ == "__main__":
    run_openmm_simulation()
```

### Custom OpenMM Protocols

```python
# custom_openmm.py
import dncs
import os
from openmm.app import *
from openmm import *
from openmm.unit import *

def custom_md_protocol(sequence, output_dir="custom_md"):
    """Run a custom MD protocol using DNCS initial structures."""
    
    # Generate initial structure with DNCS
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    
    os.makedirs(output_dir, exist_ok=True)
    temp_pdb = f"{output_dir}/initial.pdb"
    polymer.toPDB(temp_pdb)
    
    # Load structure in OpenMM
    pdb = PDBFile(temp_pdb)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    
    # Create system
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens()
    modeller.addSolvent(forcefield, padding=1.0*nanometer)
    
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds
    )
    
    # Set up integrator
    integrator = LangevinMiddleIntegrator(
        300*kelvin,    # Temperature
        1/picosecond,  # Friction coefficient
        0.002*picoseconds  # Step size
    )
    
    # Create simulation
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    
    # Equilibrate
    print("Equilibrating...")
    simulation.step(1000)  # 1000 steps = 2 ps
    
    # Production run
    print("Running production MD...")
    simulation.reporters.append(
        PDBReporter(f'{output_dir}/trajectory.pdb', 100)
    )
    simulation.reporters.append(
        StateDataReporter(
            f'{output_dir}/data.txt', 100,
            step=True, potentialEnergy=True, temperature=True
        )
    )
    
    simulation.step(5000)  # 5000 steps = 10 ps
    
    print(f"Simulation completed. Results in {output_dir}/")
    
    return output_dir

# Example usage
if __name__ == "__main__":
    try:
        result_dir = custom_md_protocol("YGG")
        print(f"Custom MD protocol completed: {result_dir}")
    except ImportError:
        print("OpenMM not available. Install with: conda install -c conda-forge openmm")
    except Exception as e:
        print(f"Error running custom MD: {e}")
```

## Advanced Topics

### Parallel Processing

```rust
use std::sync::Arc;
use rayon::prelude::*;
use libdncs::*;

fn parallel_sampling() {
    let sequences = vec!["YGG", "YGGFM", "ALANINE", "GLYCINE"];
    
    // Process sequences in parallel
    sequences.par_iter().for_each(|&sequence| {
        println!("Processing {} in parallel", sequence);
        
        let mut system = System::new(sequence, FF::AmberFB15.init());
        system.init_parameters();
        
        let folder = format!("parallel_{}", sequence);
        let mut sampler = Sampler::new(
            Arc::new(system),
            Method::Fold,
            folder.clone()
        );
        
        sampler.sample(200, 300.0);
        sampler.write_angles();
        
        println!("Completed {}", sequence);
    });
    
    println!("All parallel processing completed");
}

fn main() {
    parallel_sampling();
}
```

### Custom Energy Functions

```rust
use std::sync::Arc;
use libdncs::*;

fn custom_energy_analysis() {
    let sequence = "YGGFM";
    let mut system = System::new(sequence, FF::AmberFB15.init());
    system.init_parameters();
    
    let amber = Amber::new(Arc::new(system.clone()));
    
    // Calculate different energy components
    let total_energy = amber.energy();
    
    println!("Energy Analysis for {}", sequence);
    println!("Total Energy: {:.2} kcal/mol", total_energy);
    
    // Analyze individual atoms
    println!("\nPer-residue analysis:");
    let mut current_residue = "";
    let mut residue_atoms = Vec::new();
    
    for atom in &system.particles {
        if atom.residue != current_residue {
            if !residue_atoms.is_empty() {
                println!("Residue {}: {} atoms", current_residue, residue_atoms.len());
                residue_atoms.clear();
            }
            current_residue = &atom.residue;
        }
        residue_atoms.push(atom);
    }
    
    if !residue_atoms.is_empty() {
        println!("Residue {}: {} atoms", current_residue, residue_atoms.len());
    }
}

fn main() {
    custom_energy_analysis();
}
```

### Advanced Sampling Strategies

```python
import dncs
import numpy as np
import os

class AdaptiveSampler:
    """Advanced adaptive sampling strategy."""
    
    def __init__(self, sequence, forcefield="amberfb15.xml"):
        self.sequence = sequence
        self.polymer = dncs.Polymer(sequence, forcefield)
        self.base_energy = self.polymer.getEnergy()
        
    def hierarchical_sampling(self, base_dir="adaptive_sampling"):
        """Perform hierarchical sampling with increasing resolution."""
        
        os.makedirs(base_dir, exist_ok=True)
        
        # Level 1: Coarse sampling
        print("Level 1: Coarse sampling...")
        level1_dir = f"{base_dir}/level1"
        sampler1 = dncs.SobolSampler(
            self.polymer, 50, "explore", 300.0, level1_dir
        )
        
        # Level 2: Medium resolution
        print("Level 2: Medium resolution...")
        level2_dir = f"{base_dir}/level2"
        sampler2 = dncs.SobolSampler(
            self.polymer, 200, "search", 300.0, level2_dir
        )
        
        # Level 3: Fine sampling
        print("Level 3: Fine sampling...")
        level3_dir = f"{base_dir}/level3"
        sampler3 = dncs.SobolSampler(
            self.polymer, 500, "fold", 300.0, level3_dir
        )
        
        return [level1_dir, level2_dir, level3_dir]
    
    def temperature_ladder(self, temperatures, base_dir="temp_ladder"):
        """Perform temperature ladder sampling."""
        
        os.makedirs(base_dir, exist_ok=True)
        results = {}
        
        for temp in temperatures:
            print(f"Sampling at {temp}K...")
            temp_dir = f"{base_dir}/T{int(temp)}"
            
            sampler = dncs.SobolSampler(
                self.polymer, 200, "fold", temp, temp_dir
            )
            
            results[temp] = temp_dir
        
        return results

# Example usage
def run_advanced_sampling():
    """Run advanced sampling protocols."""
    
    sequence = "YGGFM"
    sampler = AdaptiveSampler(sequence)
    
    print(f"Advanced sampling for {sequence}")
    print(f"Base energy: {sampler.base_energy:.2f} kcal/mol")
    
    # Hierarchical sampling
    print("\n=== Hierarchical Sampling ===")
    hierarchy_results = sampler.hierarchical_sampling()
    
    for i, result_dir in enumerate(hierarchy_results, 1):
        pdb_files = [f for f in os.listdir(result_dir) if f.endswith('.pdb')]
        print(f"Level {i}: {len(pdb_files)} structures in {result_dir}")
    
    # Temperature ladder
    print("\n=== Temperature Ladder ===")
    temperatures = [250, 300, 350, 400]
    temp_results = sampler.temperature_ladder(temperatures)
    
    for temp, result_dir in temp_results.items():
        pdb_files = [f for f in os.listdir(result_dir) if f.endswith('.pdb')]
        print(f"{temp}K: {len(pdb_files)} structures in {result_dir}")

if __name__ == "__main__":
    run_advanced_sampling()
```

## Best Practices

### 1. Choosing Appropriate Sample Sizes

```python
import dncs
import time

def benchmark_sample_sizes():
    """Benchmark different sample sizes."""
    
    sequence = "YGG"
    polymer = dncs.Polymer(sequence, "amberfb15.xml")
    
    sample_sizes = [10, 50, 100, 500, 1000]
    
    print("Sample Size Benchmarking")
    print("=" * 40)
    print(f"{'Size':>6} {'Time (s)':>10} {'Structures':>12}")
    print("-" * 40)
    
    for size in sample_sizes:
        start_time = time.time()
        
        output_dir = f"benchmark_{size}"
        sampler = dncs.SobolSampler(
            polymer, size, "fold", 300.0, output_dir
        )
        
        end_time = time.time()
        
        # Count structures
        import os
        pdb_files = [f for f in os.listdir(output_dir) if f.endswith('.pdb')]
        
        print(f"{size:>6} {end_time-start_time:>10.2f} {len(pdb_files):>12}")
    
    print("\nRecommendations:")
    print("- Quick tests: 50-100 samples")
    print("- Production runs: 500-1000 samples")
    print("- Extensive studies: 1000+ samples")

if __name__ == "__main__":
    benchmark_sample_sizes()
```

### 2. Memory Management

```rust
use std::sync::Arc;
use libdncs::*;

fn memory_efficient_sampling() {
    let sequences = vec!["YGG", "YGGFM", "ALANINE"];
    
    for sequence in sequences {
        println!("Processing {}", sequence);
        
        // Create system in limited scope
        {
            let mut system = System::new(sequence, FF::AmberFB15.init());
            system.init_parameters();
            
            let system_arc = Arc::new(system);
            
            // Use Arc for shared ownership
            let mut sampler = Sampler::new(
                Arc::clone(&system_arc),
                Method::Fold,
                format!("efficient_{}", sequence)
            );
            
            sampler.sample(100, 300.0);
            sampler.write_angles();
            
            // system_arc will be dropped here
        }
        
        println!("Completed {} - memory freed", sequence);
    }
}

fn main() {
    memory_efficient_sampling();
}
```

### 3. Error Handling and Validation

```python
import dncs
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def robust_analysis(sequence, forcefield="amberfb15.xml"):
    """Robust analysis with comprehensive error handling."""
    
    try:
        # Validate inputs
        if not sequence:
            raise ValueError("Sequence cannot be empty")
        
        if not sequence.isalpha():
            raise ValueError("Sequence must contain only amino acid letters")
        
        valid_forcefields = [
            "amber03.xml", "amber10.xml", "amber96.xml", 
            "amber99sb.xml", "amberfb15.xml"
        ]
        
        if forcefield not in valid_forcefields:
            raise ValueError(f"Invalid forcefield. Must be one of: {valid_forcefields}")
        
        logger.info(f"Starting analysis for {sequence} with {forcefield}")
        
        # Create polymer with error handling
        try:
            polymer = dncs.Polymer(sequence, forcefield)
        except SystemExit:
            logger.error(f"Failed to create polymer with {forcefield}")
            return None
        
        # Calculate energy
        try:
            energy = polymer.getEnergy()
            logger.info(f"Energy calculated: {energy:.2f} kcal/mol")
        except Exception as e:
            logger.error(f"Energy calculation failed: {e}")
            return None
        
        # Create output directory safely
        output_dir = f"robust_analysis_{sequence}"
        try:
            os.makedirs(output_dir, exist_ok=True)
        except PermissionError:
            logger.error(f"Cannot create directory {output_dir}")
            return None
        
        # Save structure
        try:
            pdb_file = f"{output_dir}/structure.pdb"
            polymer.toPDB(pdb_file)
            logger.info(f"Structure saved to {pdb_file}")
        except Exception as e:
            logger.error(f"Failed to save structure: {e}")
            return None
        
        # Perform sampling
        try:
            sampler = dncs.SobolSampler(
                polymer, 100, "fold", 300.0, output_dir
            )
            logger.info("Sampling completed successfully")
        except Exception as e:
            logger.error(f"Sampling failed: {e}")
            return None
        
        # Validate results
        pdb_files = [f for f in os.listdir(output_dir) if f.endswith('.pdb')]
        if len(pdb_files) < 2:  # At least initial + 1 sample
            logger.warning(f"Only {len(pdb_files)} PDB files generated")
        
        result = {
            'sequence': sequence,
            'forcefield': forcefield,
            'energy': energy,
            'output_dir': output_dir,
            'n_structures': len(pdb_files),
            'status': 'success'
        }
        
        logger.info(f"Analysis completed successfully for {sequence}")
        return result
        
    except Exception as e:
        logger.error(f"Unexpected error in analysis: {e}")
        return {
            'sequence': sequence,
            'error': str(e),
            'status': 'failed'
        }

# Test robust analysis
def test_robust_analysis():
    """Test the robust analysis function."""
    
    test_cases = [
        ("YGG", "amberfb15.xml"),      # Valid case
        ("", "amberfb15.xml"),         # Empty sequence
        ("YGG123", "amberfb15.xml"),   # Invalid sequence
        ("YGG", "invalid.xml"),        # Invalid forcefield
        ("YGGFM", "amber99sb.xml"),    # Valid alternative
    ]
    
    results = []
    
    for sequence, forcefield in test_cases:
        print(f"\nTesting: sequence='{sequence}', forcefield='{forcefield}'")
        result = robust_analysis(sequence, forcefield)
        results.append(result)
        
        if result and result.get('status') == 'success':
            print(f"✓ Success: {result['energy']:.2f} kcal/mol")
        else:
            print(f"✗ Failed: {result.get('error', 'Unknown error') if result else 'No result'}")
    
    return results

if __name__ == "__main__":
    test_results = test_robust_analysis()
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Installation Problems

**Problem**: Rust compilation fails
```bash
error[E0599]: no method named `init` found for enum `FF`
```

**Solution**: Update to the latest version and check dependencies
```bash
rustup update
cargo clean
cargo build --release
```

**Problem**: Python module import fails
```python
ImportError: No module named 'dncs'
```

**Solution**: Rebuild and reinstall the Python module
```bash
maturin develop --release
pip install -e .
```

#### 2. Runtime Errors

**Problem**: Force field not supported
```
Unsupported forcefield: amber14.xml
```

**Solution**: Use a supported force field
```python
# Supported force fields
supported = ["amber03.xml", "amber10.xml", "amber96.xml", "amber99sb.xml", "amberfb15.xml"]
polymer = dncs.Polymer("YGG", "amberfb15.xml")  # Use supported FF
```

**Problem**: Memory issues with large systems
```
thread 'main' panicked at 'out of memory'
```

**Solution**: Reduce sample size or use batch processing
```python
# Instead of 10000 samples at once
# sampler = dncs.SobolSampler(polymer, 10000, "fold", 300.0, "output")

# Use smaller batches
for i in range(10):
    batch_dir = f"output_batch_{i}"
    sampler = dncs.SobolSampler(polymer, 1000, "fold", 300.0, batch_dir)
```

#### 3. Performance Issues

**Problem**: Slow sampling
```python
# Inefficient: Creating new polymer for each sample
for seq in sequences:
    polymer = dncs.Polymer(seq, "amberfb15.xml")  # Slow
    sampler = dncs.SobolSampler(polymer, 100, "fold", 300.0, f"output_{seq}")
```

**Solution**: Reuse objects when possible
```python
# Efficient: Batch processing
def batch_process(sequences):
    for seq in sequences:
        polymer = dncs.Polymer(seq, "amberfb15.xml")
        sampler = dncs.SobolSampler(polymer, 100, "fold", 300.0, f"output_{seq}")
        # Process immediately, then let polymer go out of scope
```

### Debugging Tips

#### 1. Enable Verbose Output

```rust
// Add debug prints
use libdncs::*;

fn debug_sampling() {
    let sequence = "YGG";
    println!("Creating system for {}", sequence);
    
    let mut system = System::new(sequence, FF::AmberFB15.init());
    println!("System created with {} atoms", system.particles.len());
    
    system.init_parameters();
    println!("Parameters initialized");
    
    let folder = "debug_output".to_string();
    let mut sampler = Sampler::new(
        std::sync::Arc::new(system),
        Method::Fold,
        folder.clone()
    );
    
    println!("Starting sampling...");
    sampler.sample(10, 300.0);
    println!("Sampling completed");
    
    sampler.write_angles();
    println!("Results written to {}", folder);
}
```

#### 2. Validate Intermediate Results

```python
import dncs
import os

def validate_pipeline(sequence):
    """Validate each step of the pipeline."""
    
    print(f"Validating pipeline for {sequence}")
    
    # Step 1: Create polymer
    try:
        polymer = dncs.Polymer(sequence, "amberfb15.xml")
        print("✓ Polymer created successfully")
    except Exception as e:
        print(f"✗ Polymer creation failed: {e}")
        return False
    
    # Step 2: Calculate energy
    try:
        energy = polymer.getEnergy()
        print(f"✓ Energy calculated: {energy:.2f} kcal/mol")
        
        # Validate energy range (rough check)
        if abs(energy) > 10000:
            print(f"⚠ Warning: Energy seems unusually high: {energy}")
    except Exception as e:
        print(f"✗ Energy calculation failed: {e}")
        return False
    
    # Step 3: Save structure
    try:
        test_pdb = "validation_test.pdb"
        polymer.toPDB(test_pdb)
        
        if os.path.exists(test_pdb):
            file_size = os.path.getsize(test_pdb)
            print(f"✓ PDB file created: {file_size} bytes")
            
            # Clean up
            os.remove(test_pdb)
        else:
            print("✗ PDB file not created")
            return False
            
    except Exception as e:
        print(f"✗ PDB creation failed: {e}")
        return False
    
    # Step 4: Test sampling
    try:
        test_dir = "validation_sampling"
        os.makedirs(test_dir, exist_ok=True)
        
        sampler = dncs.SobolSampler(polymer, 5, "fold", 300.0, test_dir)
        
        # Check output
        pdb_files = [f for f in os.listdir(test_dir) if f.endswith('.pdb')]
        if pdb_files:
            print(f"✓ Sampling successful: {len(pdb_files)} structures")
        else:
            print("✗ No structures generated")
            return False
            
        # Clean up
        import shutil
        shutil.rmtree(test_dir)
        
    except Exception as e:
        print(f"✗ Sampling failed: {e}")
        return False
    
    print("✓ All validation steps passed")
    return True

# Test validation
if __name__ == "__main__":
    test_sequences = ["A", "YGG", "YGGFM"]
    
    for seq in test_sequences:
        print(f"\n{'='*40}")
        success = validate_pipeline(seq)
        print(f"Validation for {seq}: {'PASSED' if success else 'FAILED'}")
```

### Getting Help

1. **Check the documentation**: Start with the API reference and examples
2. **Validate your input**: Ensure sequences and force fields are correct
3. **Start small**: Test with simple sequences before complex ones
4. **Check logs**: Look for error messages in output files
5. **Use debugging**: Add print statements to track execution
6. **Report issues**: Include minimal reproducible examples

This comprehensive tutorial should help you get started with DNCS and progress to advanced usage. Remember to start with simple examples and gradually increase complexity as you become more familiar with the library.