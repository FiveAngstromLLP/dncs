# Digital Nets Conformational Sampling (DNCS)

DNCS is a tool for molecular conformational sampling using digital nets. It provides both a command-line interface and Python bindings for generating and analyzing molecular conformations.

## Installation

### Prerequisites

- Rust toolchain (1.70.0 or later)
- Cargo (Rust's package manager)
- Python 3.7+ (for Python bindings)
- Maturin (for building Python bindings)

### Installing Rust

**For Unix-based systems (Linux, macOS):**
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**For Windows:**
- Download and run [rustup-init.exe](https://rustup.rs)


### Installing from Source

1. Clone the repository:
```bash
git clone https://github.com/FiveAngstromLLP/dncs.git
cd dncs
```

2. Install the command-line tool:
```bash
cargo install --path .
```

3. (Optional) Install Python bindings:
```bash
pip install maturin
cd python
maturin develop --release
```

## Usage

### Command-line Interface

1. Generate a configuration file:
```bash
dncs -c
```
This creates a `dncs.json` file with some example parameters:
```json
{
  "Generate": {
    "molecule": "Sample",
    "sequence": "YGGFM",
    "n_samples": 10,
    "forcefield": "amberfb15.xml",
    "minimize": true,
    "grid": 4
  }
}
```

2. Run simulation with command-line arguments:
```bash
dncs -N molecule_name -s SEQUENCE -n NUMBER_OF_SAMPLES -f FORCEFIELD [-m] [-g GRID_SIZE]
```

Example:
```bash
dncs -N test -s YGGFM -n 100 -f amberfb15.xml -m -g 4
```

3. Run simulation using config file:
```bash
dncs
```

### Command-line Options

- `-N, --molecule`: Molecule name
- `-s, --sequence`: Amino acid sequence
- `-n, --samples`: Number of samples to generate
- `-f, --forcefield`: Force field selection
  - Available options:
    - amber03.xml
    - amber10.xml
    - amber96.xml
    - amber99sb.xml
    - amberfb15.xml
- `-m, --minimize`: Enable energy minimization
- `-g, --grid`: Grid size for sample division
- `-c, --config`: Generate default configuration file

### Python Interface

```python
import dncs

# Create a polymer system
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# Get energy
energy = polymer.getEnergy()

# Generate PDB file
polymer.toPDB("output.pdb")

# Generate conformational samples
sampler = dncs.SobolSampler(polymer, n_samples=100, grid=4)

# Write output files
sampler.write_angles("angles.out")
sampler.toPDB("conformations.pdb")
sampler.toPDBFiles("conformations")
```

## Output Files

The program creates a `Result` directory with the following structure:
```
Result/
└── molecule_name/
    ├── Sampled.out       # Initial sampled angles
    ├── Sampled.pdb       # Initial sampled structures
    ├── molecule_name.out # Minimized angles (if -m is used)
    └── molecule_name.pdb # Minimized structures (if -m is used)
```
