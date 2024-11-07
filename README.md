# DNCS 1.0 (Digital Nets Conformational Sampling)

DNCS 1.0 is an Enhanced Conformational Sampling tool using digital nets. It provides both a command-line interface and Python bindings for generating and analyzing molecular conformations.The DNCS tool is based on the following publications.
We request you to cite us when using our tool.

```bibtex
@Article{D4CP01891E,
author ="J., Abraham Rebairo and D., Sam Paul and Arumainathan, Stephen",
title  ="Digital nets conformational sampling (DNCS) – an enhanced sampling technique to explore the conformational space of intrinsically disordered peptides",
journal  ="Phys. Chem. Chem. Phys.",
year  ="2024",
volume  ="26",
issue  ="34",
pages  ="22640-22655",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/D4CP01891E",
url  ="http://dx.doi.org/10.1039/D4CP01891E",
abstract  ="We propose digital nets conformational sampling (DNCS) – an enhanced sampling technique to explore the conformational ensembles of peptides{,} especially intrinsically disordered peptides (IDPs). The DNCS algorithm relies on generating history-dependent samples of dihedral variables using bitwise XOR operations and binary angle measurements (BAM). The algorithm was initially studied using met-enkephalin{,} a highly elusive neuropeptide. The DNCS method predicted near-native structures and the energy landscape of met-enkephalin was observed to be in direct correlation with earlier studies on the neuropeptide. Clustering analysis revealed that there are only 24 low-lying conformations of the molecule. The DNCS method has then been tested for predicting optimal conformations of 42 oligopeptides of length varying from 3 to 8 residues. The closest-to-native structures of 86% of cases are near-native and 24% of them have a root mean square deviation of less than 1.00 Å with respect to their crystal structures. The results obtained reveal that the DNCS method performs well{,} that too in less computational time."}


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
    "molecule": "Met-Enkephalin",
    "sequence": "YGGFM",
    "n_samples": 10,
    "forcefield": "amberfb15.xml",
    "minimize": true,
    "grid": 4
  }
}
```

2. Run Sampling with command-line arguments:
```bash
dncs -N molecule_name -s SEQUENCE -n NUMBER_OF_SAMPLES -f FORCEFIELD [-m] [-g GRID_SIZE]
```

Example:
```bash
dncs -N test -s YGGFM -n 100 -f amberfb15.xml -m -g 4
```

3. Run Sampling using config file:
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

## Integration with openmm

You can find a template for OpenMM integration in the https://github.com/FiveAngstromLLP/dncs-Sampling.git repository.
