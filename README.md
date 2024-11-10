[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2FD4CP01891E-blue)](https://doi.org/10.1039/D4CP01891E)
[![DOI](https://zenodo.org/badge/863510320.svg)](https://doi.org/10.5281/zenodo.14054733)
# DNCS 1.0 (Digital Nets Conformational Sampling)

DNCS 1.0 is an Enhanced Conformational Sampling tool using digital nets. It provides  Python bindings for enhanced sampling of molecular conformations.


## Prerequisites

Before installing DNCS, ensure you have:

- Rust toolchain (1.70.0 or later)
- Cargo (Rust's package manager)
- Python 3.7+ (for Python bindings)
- Maturin (for building Python bindings)

## Installation

**For Unix-based systems (Linux, macOS):**
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**For Windows:**
- Download and run [rustup-init.exe](https://rustup.rs)

### Usage

- Installation:

  ```bash
  just install
  ```
  - Installs Python requirements
  - Installs maturin
  - Builds Rust components

- Running:
  ```bash
  just run
  ```
  - Executes main DNCS script

Note: Uses DNCS_FOLDER environment variable for path management

## Configuration

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
    ├── dncs.log                # Log file
    ├── Langevin/               # Langevin dynamics results
    │   ├── Equilibrated_*.pdb  # Equilibrated structures
    │   └── equilibrated.out    # Angle measurements
    ├── MDSimulation/           # Molecular dynamics results
    │   └── simulated_*.pdb     # Simulation trajectory structures
    ├── Minimized/              # Energy minimization results
    │   ├── Minimized_*.pdb     # Minimized structures
    │   └── minimized.out       # Final angles
    ├── Sampled/                # Initial sampling results
    │   ├── angles.out          # Initial angles
    │   └── sample_*.pdb        # Initial structures
    ├── sampled.pdb             # Sampled structure
    ├── minimized.pdb           # Minimized structure
    ├── equilibrated.pdb        # Equilibrated structure
    └── linear.pdb              # Linear structure
```
