[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2FD4CP01891E-blue)](https://doi.org/10.1039/D4CP01891E)

# DNCS 1.0 (Digital Nets Conformational Sampling)

DNCS 1.0 is an Enhanced Conformational Sampling tool using digital nets. It provides both a command-line interface and Python bindings for generating and analyzing molecular conformations.


## Prerequisites

Before installing DNCS, ensure you have:

- Rust toolchain (1.70.0 or later)
- Cargo (Rust's package manager)
- Python 3.7+ (for Python bindings)
- Maturin (for building Python bindings)

## Installation

### 1. Installing Rust

**For Unix-based systems (Linux, macOS):**
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**For Windows:**
- Download and run [rustup-init.exe](https://rustup.rs)

### 2. Installing DNCS

1. Clone the repository:
```bash
git clone https://github.com/FiveAngstromLLP/dncs.git
cd dncs
```

2. Install the command-line tool:
```bash
cargo install --path .
```

3. Install Python bindings (optional):
```bash
pip install maturin
cd python
maturin develop --release
```

## Configuration

### Rust CLI Configuration Options

#### JSON Configuration File

Create a configuration file for CLI sampling:
```bash
dncs -c
```

This creates `dncs.json` with parameters for enhanced conformational sampling:
```json
{
  "Generate": {
    "molecule": "Met-Enkephalin",
    "sequence": "YGGFM",
    "n_samples": 10,
    "forcefield": "amberfb15.xml",
    "minimize": true,
    "grid": 4
  }
}
```

#### Command-line Options

Basic usage:
```bash
dncs -N molecule_name -s SEQUENCE -n NUMBER_OF_SAMPLES -f FORCEFIELD [-m] [-g GRID_SIZE]
```

Example:
```bash
dncs -N test -s YGGFM -n 100 -f amberfb15.xml -m -g 4
```

Run with config file:
```bash
dncs
```

CLI Options:
- `-N, --molecule`: Molecule name
- `-s, --sequence`: Amino acid sequence
- `-n, --samples`: Number of samples to generate
- `-f, --forcefield`: Force field selection (amber03.xml|amber10.xml|amber96.xml|amber99sb.xml|amberfb15.xml)
- `-m, --minimize`: Enable energy minimization
- `-g, --grid`: Grid size for adaptive sampling.
- `-c, --config`: Generate default configuration file

### Output Structure [Rust Sampling]

Results are stored in:
```
Result/
└── molecule_name/
    ├── Sampled.out       # Initial sampled angles
    ├── Sampled.pdb       # Initial sampled structures
    ├── molecule_name.out # Minimized angles (if -m is used)
    └── molecule_name.pdb # Minimized structures (if -m is used)
```

### Advanced Configuration for OpenMM Simulation

For Python OpenMM simulation, create a `dncs.toml` file:

```toml
[simulation]
moleculename = "Met-Enkephalin" #MOLECULE NAME
sequence = "YGGFM" #INPUT SEQUENCE
interface = "openmm" #INTERFACE TO OPENMM (OpenMM Should be installed from https://openmm.org/
n_samples = 100 # Number of Samples (This divides the timesteps into segments)
temp = 300.0 # NVT equilibration temperature.
forcefield = ["amber14.xml", "amber14/tip3pfb.xml"] #Force Field specification
device = "CUDA" #Device to run MD simulation
solvent = 10000 # Solvation
steps = 5000 # Equlibration timesteps (uses Langevin Integrator from OpenMM)
gamma = 1.0 # Friction coefficient 
dt = 0.002 # Integrator timestep
md_steps = 5000 # MD timestep
grid = 4 # Adaptive Sampling
```

### Output Structure [Python Simulation]

The Python interface generates the following directory structure:
```
Result/moleculename/
    ├── dncs.log                    # Log file
    ├── Langevin/                   # Langevin dynamics results
    │   ├── Equilibrated_*.pdb      # Equilibrated structures
    │   └── equilibrated.out        # Angle measurements
    ├── MDSimulation/               # Molecular dynamics results
    │   └── simulated_*.pdb         # Simulation trajectory structures
    ├── Minimized/                  # Energy minimization results
    │   ├── Minimized_*.pdb        # Minimized structures
    │   └── minimized.out          # Final angles
    ├── Sampled/                    # Initial sampling results
    │   ├── angles.out             # Initial angles
    │   └── sample_*.pdb           # Initial structures
    ├── multi_structure.pdb         # Combined structure file
    ├── result.pdb                 # Final structure
    └── structure.pdb              # Input structure
```

## Usage

### Using the Justfile

For development tasks:

Installation:
```bash
just install
```
- Installs Python requirements
- Installs maturin
- Builds Rust components

Running:
```bash
just run
```
- Executes main DNCS script

Note: Uses DNCS_FOLDER environment variable for path management
