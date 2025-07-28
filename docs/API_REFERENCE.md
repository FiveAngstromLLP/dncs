# DNCS API Reference

## Table of Contents

1. [Overview](#overview)
2. [Rust API](#rust-api)
   - [Core Structs](#core-structs)
   - [Force Fields](#force-fields)
   - [Sampling](#sampling)
   - [Minimization](#minimization)
   - [Parsing](#parsing)
3. [Python API](#python-api)
   - [Functions](#functions)
   - [Classes](#classes)
4. [Command Line Interface](#command-line-interface)
5. [Configuration](#configuration)
6. [Examples](#examples)

## Overview

DNCS (Digital Nets Conformational Sampling) is a molecular simulation library that provides enhanced conformational sampling using digital nets. It offers both Rust and Python APIs for molecular system manipulation, energy calculations, and conformational sampling.

## Rust API

### Core Structs

#### `System`

The primary struct representing a molecular system with particles and force field parameters.

```rust
pub struct System {
    pub seq: String,
    pub forcefield: ForceField,
    pub particles: Particles,
    pub dihedral: Vec<(Atom, Atom, Atom, Atom)>,
    pub dihedral_angle: Vec<(Atom, Atom, Atom, Atom)>,
    pub firstbonded: Vec<Particles>,
    pub secondbonded: Vec<Particles>,
    pub nonbonded: Vec<Particles>,
    pub bonded1_4: Vec<Particles>,
    pub hydrogen: Vec<(Atom, Atom)>,
}
```

**Methods:**

- `new(seq: &str, forcefield: ForceField) -> Self`
  
  Creates a new System from an amino acid sequence string.
  
  ```rust
  use libdncs::*;
  
  let system = System::new("YGGFM", FF::AmberFB15.init());
  ```

- `from_pdb(file: &str, forcefield: ForceField) -> Self`
  
  Creates a System by parsing a PDB file.
  
  ```rust
  let system = System::from_pdb("protein.pdb", FF::AmberFB15.init());
  ```

- `init_parameters(&mut self)`
  
  Initializes all system parameters including atom types, neighbors, and energy parameters.
  
  ```rust
  let mut system = System::new("YGGFM", FF::AmberFB15.init());
  system.init_parameters();
  ```

- `get_atomtype(&mut self)`
  
  Assigns atom types based on force field definitions.

- `get_dihedral(&mut self)`
  
  Calculates dihedral angles in the system.

- `get_neighbours(&mut self)`
  
  Determines bonded and non-bonded neighbors for each atom.

- `to_pdb(&self, filename: &str)`
  
  Exports the system to a PDB file.
  
  ```rust
  system.to_pdb("output.pdb");
  ```

- `export_pdb(&self) -> String`
  
  Returns the system as a PDB format string.

- `dihedral_log(&self, foldername: &str)`
  
  Logs dihedral angles to a specified folder.

#### `Atom`

Represents an individual atom with all its properties.

```rust
pub struct Atom {
    pub record: String,
    pub serial: usize,
    pub name: String,
    pub residue: String,
    pub chain: String,
    pub residue_number: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub temperature_factor: f64,
    pub element: String,
    pub charge: Option<f64>,
    pub epsilon: Option<f64>,
    pub sigma: Option<f64>,
    pub typeid: Option<String>,
    pub atomtype: Option<String>,
}
```

**Methods:**

- `new(line: String) -> Self`
  
  Creates an Atom from a PDB line string.

### Force Fields

#### `FF` Enum

Represents available force fields.

```rust
pub enum FF {
    Amber03,
    Amber10,
    Amber96,
    Amber99SB,
    AmberFB15,
}
```

**Methods:**

- `init(&self) -> ForceField`
  
  Initializes the force field parameters.
  
  ```rust
  let ff = FF::AmberFB15.init();
  ```

#### `Amber`

Force field implementation for energy calculations.

```rust
pub struct Amber {
    pub system: Arc<System>,
}
```

**Methods:**

- `new(system: Arc<System>) -> Self`
  
  Creates a new Amber force field instance.
  
  ```rust
  use std::sync::Arc;
  
  let amber = Amber::new(Arc::new(system));
  ```

- `energy(&self) -> f64`
  
  Calculates the total energy of the system in kcal/mol.
  
  ```rust
  let total_energy = amber.energy();
  println!("Energy: {} kcal/mol", total_energy);
  ```

### Sampling

#### `Sobol`

Generates Sobol sequences for uniform sampling in high-dimensional spaces.

```rust
pub struct Sobol {
    count: usize,
    total: usize,
    current: Vec<usize>,
    direction: Vec<Vec<usize>>,
}
```

**Methods:**

- `new(dimension: usize) -> Self`
  
  Creates a new Sobol sequence generator.
  
  ```rust
  let sobol = Sobol::new(10); // 10-dimensional space
  for point in sobol.take(100) {
      println!("{:?}", point); // Vec<f64> with values in [0,1]
  }
  ```

#### `RotateAtDihedral`

Handles rotation of molecular systems around dihedral angles.

```rust
pub struct RotateAtDihedral {
    system: System,
    rotated: Particles,
}
```

**Methods:**

- `new(system: Arc<System>) -> Self`
  
  Creates a new rotation handler.
  
  ```rust
  let rotator = RotateAtDihedral::new(Arc::new(system));
  ```

- `from_pdb(file: &str) -> Vec<f64>`
  
  Extracts dihedral angles from a PDB file.
  
  ```rust
  let angles = RotateAtDihedral::from_pdb("structure.pdb");
  ```

- `rotate(&mut self, angle: Vec<f64>)`
  
  Rotates the system according to specified dihedral angles.
  
  ```rust
  let angles = vec![1.0, 2.0, 3.0]; // radians
  rotator.rotate(angles);
  ```

- `rotated_energy(&mut self, angle: Vec<f64>) -> f64`
  
  Calculates energy after rotation.

- `current_dihedral(&self) -> Vec<f64>`
  
  Returns current dihedral angles.

- `energy(&self) -> f64`
  
  Returns energy of current configuration.

- `to_pdbstring(&self, model: usize, energy: f64) -> String`
  
  Converts to PDB format string.

- `dihedral_angle(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> f64`
  
  Calculates dihedral angle between four atoms.

#### `Method` Enum

Sampling methods available.

```rust
pub enum Method {
    Fold,
    Search,
    Explore,
}
```

#### `Sampler`

Main sampling class for conformational sampling.

```rust
pub struct Sampler {
    system: Arc<System>,
    rotate: RotateAtDihedral,
    method: Method,
    folder: String,
    angles: Vec<Vec<f64>>,
    sample: Vec<System>,
    energy: Vec<f64>,
}
```

**Methods:**

- `new(system: Arc<System>, method: Method, folder: String) -> Self`
  
  Creates a new sampler.
  
  ```rust
  let sampler = Sampler::new(
      Arc::new(system),
      Method::Fold,
      "output_folder".to_string()
  );
  ```

- `sample(&mut self, max: usize, temp: f64)`
  
  Generates conformational samples.
  
  ```rust
  sampler.sample(1000, 300.0); // 1000 samples at 300K
  ```

- `write_angles(&self)`
  
  Writes sampled angles to file.

- `transform_angle(&self, angle: Vec<f64>) -> Vec<f64>`
  
  Transforms angles based on sampling method.

### Minimization

#### `Minimizer`

Energy minimization using L-BFGS algorithm.

```rust
pub struct Minimizer {
    pub samplefolder: String,
    pub folder: String,
    pub top_n_sample: usize,
    pub forcefield: ForceField,
    pub sample: Option<RotateAtDihedral>,
    pub energy: Vec<(usize, f64)>,
    pub angles: Vec<Vec<f64>>,
}
```

**Methods:**

- `new(folder: String, forcefield: ForceField, top_n_sample: usize) -> Self`
  
  Creates a new minimizer.
  
  ```rust
  let minimizer = Minimizer::new(
      "output".to_string(),
      FF::AmberFB15.init(),
      10
  );
  ```

- `minimize_all(&mut self)`
  
  Minimizes all sampled structures.

- `minimize(&mut self, system: Arc<System>, model: usize)`
  
  Minimizes a specific system.

- `conformational_sort(&mut self, temp: f64)`
  
  Sorts conformations by Boltzmann weights.

- `write_angles(&self)`
  
  Writes minimized angles to file.

- `rename(&self)`
  
  Renames output files based on energy ranking.

### Parsing

#### Parser Functions

- `generate(seq: &str) -> Vec<Atom>`
  
  Generates atoms from amino acid sequence.
  
  ```rust
  let atoms = parser::generate("YGGFM");
  ```

- `atoms_to_pdbstring(atoms: Vec<Atom>) -> String`
  
  Converts atoms to PDB format string.

- `pdb_to_atoms(pdb_string: &str) -> Vec<Atom>`
  
  Parses PDB string to atoms.

- `atoms_to_seq(atoms: Vec<Atom>) -> String`
  
  Converts atoms back to sequence string.

## Python API

### Functions

#### `getPDB(seq: str, filename: str) -> None`

Generates a PDB file from an amino acid sequence.

```python
import dncs

dncs.getPDB("YGGFM", "output.pdb")
```

#### `pdb_to_angle(filename: str) -> str`

Extracts dihedral angles from a PDB file as a comma-separated string.

```python
angles_csv = dncs.pdb_to_angle("structure.pdb")
print(angles_csv)  # "1.23, 2.45, 3.67, ..."
```

### Classes

#### `Polymer`

Represents a polymer system with force field parameters.

```python
class Polymer:
    def __init__(self, seq: str, forcefield: str) -> None
    def getEnergy(self) -> float
    def toPDB(self, filename: str) -> None
    def dihedral(self, foldername: str) -> None
```

**Constructor:**

```python
# Supported forcefields: amber03.xml, amber10.xml, amber96.xml, amber99sb.xml, amberfb15.xml
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")
```

**Methods:**

- `getEnergy() -> float`
  
  Calculates system energy in kcal/mol.
  
  ```python
  energy = polymer.getEnergy()
  print(f"Energy: {energy} kcal/mol")
  ```

- `toPDB(filename: str) -> None`
  
  Saves structure to PDB file.
  
  ```python
  polymer.toPDB("structure.pdb")
  ```

- `dihedral(foldername: str) -> None`
  
  Logs dihedral angles to folder.
  
  ```python
  polymer.dihedral("angles_output")
  ```

#### `SobolSampler`

Sobol sequence sampler for conformational sampling.

```python
class SobolSampler:
    def __init__(self, system: Polymer, no_of_samples: int, method: str, temp: float, folder: str) -> None
```

**Constructor:**

```python
sampler = dncs.SobolSampler(
    polymer,           # Polymer system
    1000,             # Number of samples
    "fold",           # Method: "fold", "search", or "explore"  
    300.0,            # Temperature in Kelvin
    "output_folder"   # Output directory
)
```

## Command Line Interface

The DNCS binary provides command-line access to sampling functionality.

```bash
# Basic usage
dncs --sequence "YGGFM" --samples 1000 --output results/

# With specific force field
dncs --sequence "YGGFM" --forcefield amber99sb --samples 500

# From PDB file
dncs --pdb input.pdb --samples 1000 --minimize
```

**Options:**
- `--sequence, -s`: Amino acid sequence
- `--pdb`: Input PDB file
- `--samples, -n`: Number of samples (default: 100)
- `--forcefield, -f`: Force field (default: amberfb15)
- `--output, -o`: Output directory
- `--temperature, -t`: Temperature in Kelvin (default: 300.0)
- `--method, -m`: Sampling method (fold/search/explore)
- `--minimize`: Enable energy minimization

## Configuration

### TOML Configuration

Create a `dncs.toml` file for Python simulations:

```toml
[simulation]
moleculename = "6RRO"
folder = "Result"
sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
interface = "openmm"
n_samples = 100
md_simulation = 10
temp = 300.0
forcefield = ["amber14-all.xml", "amber14/tip3p.xml"]
device = "HIP"
solvent = 10
steps = 50
gamma = 1.0
dt = 0.002
md_steps = 1
method = "fold"
```

**Parameters:**
- `moleculename`: Molecule identifier
- `sequence`: Amino acid sequence
- `n_samples`: Number of conformational samples
- `md_simulation`: Number of top samples for MD simulation
- `temp`: Temperature in Kelvin
- `forcefield`: OpenMM force field files
- `device`: Computation device (CPU/CUDA/HIP)
- `solvent`: Number of solvent molecules
- `method`: Sampling method

## Examples

### Complete Rust Example

```rust
use std::sync::Arc;
use libdncs::*;

fn main() {
    // 1. Create a molecular system
    let mut system = System::new("YGGFM", FF::AmberFB15.init());
    system.init_parameters();
    
    // 2. Calculate initial energy
    let amber = Amber::new(Arc::new(system.clone()));
    let initial_energy = amber.energy();
    println!("Initial energy: {} kcal/mol", initial_energy);
    
    // 3. Perform conformational sampling
    let mut sampler = Sampler::new(
        Arc::new(system),
        Method::Fold,
        "results".to_string()
    );
    sampler.sample(1000, 300.0);
    sampler.write_angles();
    
    // 4. Energy minimization
    let mut minimizer = Minimizer::new(
        "results".to_string(),
        FF::AmberFB15.init(),
        10
    );
    minimizer.minimize_all();
    minimizer.conformational_sort(300.0);
    minimizer.write_angles();
}
```

### Complete Python Example

```python
import dncs
from integrator import DncsIntegrator
import toml

# 1. Create polymer system
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# 2. Calculate energy
energy = polymer.getEnergy()
print(f"Initial energy: {energy} kcal/mol")

# 3. Generate conformational samples
sampler = dncs.SobolSampler(
    polymer,
    1000,      # samples
    "fold",    # method
    300.0,     # temperature
    "output"   # folder
)

# 4. Save structure
polymer.toPDB("initial_structure.pdb")
polymer.dihedral("dihedral_angles")

# 5. Run MD simulation (requires OpenMM)
config = toml.load("dncs.toml")["simulation"]
integrator = DncsIntegrator(config)
integrator.run_integrator()
```

### Sobol Sampling Example

```rust
use libdncs::*;

fn main() {
    // Generate 100 points in 5D space
    let sobol = Sobol::new(5);
    let points: Vec<Vec<f64>> = sobol.take(100).collect();
    
    for (i, point) in points.iter().enumerate() {
        println!("Point {}: {:?}", i, point);
    }
}
```

### Dihedral Rotation Example

```rust
use std::sync::Arc;
use libdncs::*;

fn main() {
    let mut system = System::new("YGGFM", FF::AmberFB15.init());
    system.init_parameters();
    
    let mut rotator = RotateAtDihedral::new(Arc::new(system));
    
    // Get current dihedral angles
    let current_angles = rotator.current_dihedral();
    println!("Current angles: {:?}", current_angles);
    
    // Rotate by 30 degrees (π/6 radians) for each dihedral
    let new_angles: Vec<f64> = current_angles
        .iter()
        .map(|_| std::f64::consts::PI / 6.0)
        .collect();
    
    rotator.rotate(new_angles);
    let new_energy = rotator.energy();
    println!("Energy after rotation: {} kcal/mol", new_energy);
}
```

## Error Handling

### Rust

Most functions return `Result` types or panic on invalid input:

```rust
// Force field validation
let ff = match FF::from_str("amber99sb.xml") {
    Ok(forcefield) => forcefield,
    Err(e) => {
        eprintln!("Invalid forcefield: {}", e);
        return;
    }
};
```

### Python

Python functions raise appropriate exceptions:

```python
try:
    polymer = dncs.Polymer("INVALID", "invalid_ff.xml")
except ValueError as e:
    print(f"Error: {e}")
```

## Performance Considerations

1. **Parallel Processing**: The library uses Rayon for parallel computations
2. **Memory Management**: Use `Arc<System>` for shared ownership
3. **Large Systems**: Consider chunking large conformational spaces
4. **Energy Calculations**: Amber force field calculations are computationally intensive

## Output Structure

DNCS generates organized output directories:

```
Result/moleculename/
├── dncs.log                # Simulation log
├── Sampled/                # Initial samples
│   ├── angles.out          # Dihedral angles
│   └── sample_*.pdb        # Sample structures
├── Minimized/              # Energy minimized
│   ├── Minimized_*.pdb     # Minimized structures
│   └── minimized.out       # Final angles
├── Langevin/               # MD equilibration
│   ├── Equilibrated_*.pdb  # Equilibrated structures
│   └── equilibrated.out    # Equilibration data
└── MDSimulation/           # Full MD simulation
    └── simulated_*.pdb     # Trajectory frames
```