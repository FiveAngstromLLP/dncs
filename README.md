<div align="center">

# 🧬 DNCS
### **Digital Nets Conformational Sampling**

*A High-Performance Peptide Structure Prediction Tool Using Quasi-Monte Carlo Methods*

[![Version](https://img. shields.io/badge/version-1. 1.2-blue.svg)](https://github.com/FiveAngstromLLP/dncs)
[![License](https://img.shields. io/badge/license-GPL%20v3-green. svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang. org/)
[![Python](https://img. shields.io/badge/python-3. 7%2B-yellow.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2FD4CP01891E-red.svg)](https://doi.org/10.1039/D4CP01891E)

---

**[Documentation](docs/TUTORIAL.md)** • **[API Reference](docs/API_REFERENCE.md)** • **[Python API](docs/PYTHON_API.md)** • **[Examples](examples/)**

</div>

---

## 🌟 Overview

**DNCS** (Digital Nets Conformational Sampling) is a cutting-edge computational tool for **peptide 3D structure prediction** that leverages the mathematical elegance of **Sobol quasi-random sequences** to achieve superior conformational space exploration compared to traditional Monte Carlo methods.

Unlike conventional random sampling approaches, DNCS employs **low-discrepancy sequences** that ensure uniform coverage of the high-dimensional dihedral angle space, enabling faster convergence to native-like peptide structures with fewer computational samples. 

<div align="center">

```
    Amino Acid Sequence                    3D Peptide Structure
         ┌─────┐                              ┌─────────────┐
         │YGGFM│  ──────────────────────────► │   🧬 PDB    │
         └─────┘        DNCS Engine           └─────────────┘
                    ┌───────────────┐
                    │ Sobol Sampling│
                    │ AMBER Energy  │
                    │ L-BFGS Optim. │
                    └───────────────┘
```

</div>

---

## ✨ Key Features

<table>
<tr>
<td width="50%">

### 🎯 **Enhanced Sampling**
- **Sobol Sequences**: Low-discrepancy quasi-Monte Carlo sampling for uniform conformational space exploration
- **Adaptive Methods**: Three specialized strategies—`Fold`, `Search`, and `Explore`
- **Up to 21,201 dimensions** supported for complex peptide systems

</td>
<td width="50%">

### ⚡ **High Performance**
- **Rust Core**: Native performance with zero-cost abstractions
- **Parallel Processing**: Multi-threaded sampling via Rayon
- **Memory Efficient**: Optimized data structures for large peptides

</td>
</tr>
<tr>
<td width="50%">

### 🔬 **Forcefields and Standard Techniques**
- **AMBER Force Fields**: FB15, 99SB, 03, 10, and 96 variants
- **Complete Energy Model**: Lennard-Jones, electrostatics, bonds, angles, and torsions
- **L-BFGS Minimization**: Gradient-based structure refinement

</td>
<td width="50%">

### 🐍 **Seamless Integration**
- **Python Bindings**: Native PyO3-based interface
- **OpenMM Compatible**: Molecular dynamics refinement
- **Standard Formats**: PDB input/output support

</td>
</tr>
</table>

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/FiveAngstromLLP/dncs.git
cd dncs

# Quick install (recommended)
just install

# Or manual installation
pip install -r requirements.txt
pip install maturin
maturin develop --release
```

### Prerequisites

| Requirement | Version | Purpose |
|-------------|---------|---------|
| **Rust** | ≥ 1. 70. 0 | Core engine compilation |
| **Python** | ≥ 3.7 | Python bindings |
| **Maturin** | Latest | Build system for PyO3 |
| **OpenMM** | Optional | MD refinement |

---

## 📖 Usage Examples

### 🐍 Python Interface

```python
import dncs

# Step 1: Create peptide system from amino acid sequence
peptide = dncs. Polymer("YGGFM", "amberfb15. xml")  # Met-enkephalin

# Step 2: Calculate initial energy
energy = peptide.getEnergy()
print(f"Initial Energy: {energy:.2f} kJ/mol")

# Step 3: Predict structures using Sobol sampling
sampler = dncs. SobolSampler(
    peptide,          # Peptide system
    1000,             # Number of conformations to generate
    "fold",           # Sampling method
    "output/"         # Results directory
)

# Step 4: Export best structure
peptide. toPDB("predicted_structure.pdb")
```

### 🦀 Rust Interface

```rust
use std::sync::Arc;
use libdncs::*;

fn main() {
    // Create peptide from sequence with AMBER FB15 force field
    let mut system = System::new("YGGFM", FF::AmberFB15. init());
    system.init_parameters();
    
    // Calculate conformational energy
    let amber = Amber::new(Arc::new(system. clone()));
    let energy = amber.energy();
    println!("Energy: {:.2} kJ/mol", energy);
    
    // Generate structure predictions
    let mut sampler = Sampler::new(
        Arc::new(system),
        Method::Fold,
        "results". to_string()
    );
    sampler.sample(1000);
}
```

---

## ⚙️ Configuration

Create a `dncs.toml` file for advanced simulations with OpenMM integration:

```toml
[simulation]
# Peptide Information
moleculename = "enkephalin"
sequence = "YGGFM"
folder = "Results"

# Sampling Configuration
method = "fold"           # Options: fold, search, explore
n_samples = 100           # Number of conformational samples
temp = 300.0              # Temperature (Kelvin)

# Force Field
forcefield = [
    "amber14-all. xml",
    "amber14/tip3p.xml",
]

# OpenMM Integration (optional)
interface = "openmm"
md_simulation = 10        # Top-N structures for MD refinement
device = "CUDA"           # CUDA, OpenCL, HIP, or CPU

# Simulation Parameters
solvent = 10              # Solvation box (Å)
steps = 50                # Equilibration steps
gamma = 1.0               # Friction coefficient (1/ps)
dt = 0.002                # Timestep (ps)
md_steps = 1              # Production MD steps
```

---

## 🎯 Sampling Methods

DNCS offers three specialized sampling strategies optimized for different prediction scenarios:

| Method | Description | Best For |
|--------|-------------|----------|
| **`fold`** | Focused backbone sampling with Ramachandran-aware angle distributions | Initial structure prediction, small peptides |
| **`search`** | Targeted conformational search around energy minima | Refinement, finding specific conformations |
| **`explore`** | Broad exploration across full conformational space | Large peptides, ensemble generation |

```python
# Example: Comparing sampling methods
for method in ["fold", "search", "explore"]:
    sampler = dncs.SobolSampler(peptide, 500, method, f"output_{method}/")
```

---

## 🔬 Force Fields

DNCS implements the complete **AMBER** force field family with accurate parameterization:

| Force Field | XML File | Description |
|-------------|----------|-------------|
| **AMBER FB15** | `amberfb15.xml` | Force Balance optimized (recommended) |
| **AMBER 99SB** | `amber99sb.xml` | Improved backbone parameters |
| **AMBER 14** | `amber14-all.xml` | Latest generation with extensive validation |
| **AMBER 03** | `amber03.xml` | Refined charge model |
| **AMBER 10** | `amber10.xml` | Updated van der Waals parameters |
| **AMBER 96** | `amber96.xml` | Original AMBER implementation |

### Energy Components

The total conformational energy is computed as:

```
E_total = E_LJ + E_elec + E_bond + E_angle + E_torsion
```

Where:
- **E_LJ**: Lennard-Jones van der Waals interactions
- **E_elec**: Coulombic electrostatic interactions  
- **E_bond**: Harmonic bond stretching
- **E_angle**: Harmonic angle bending
- **E_torsion**: Periodic dihedral torsions

---

## 📁 Output Structure

After running a prediction, DNCS generates organized results:

```
Results/moleculename/
├── 📄 dncs. log                 # Detailed simulation log
├── 📄 linear. pdb               # Initial extended structure
├── 📄 sampled.pdb              # Best sampled conformation
├── 📄 minimized.pdb            # Energy-minimized structure
├── 📄 equilibrated.pdb         # Thermally equilibrated structure
│
├── 📂 Sampled/                 # Raw sampling results
│   ├── sample_*. pdb            # Individual conformations
│   └── angles. out              # Dihedral angle data
│
├── 📂 Minimized/               # Energy minimization
│   ├── Minimized_*.pdb         # Refined structures
│   └── minimized.out           # Final angles
│
├── 📂 Langevin/                # Thermal equilibration
│   ├── Equilibrated_*.pdb      # Equilibrated conformations
│   └── equilibrated.out        # Post-equilibration data
│
└── 📂 MDSimulation/            # MD refinement (if enabled)
    └── simulated_*.pdb         # Production trajectories
```

---

## 🏗️ Architecture

```
dncs/
├── src/                        # Rust core library
│   ├── system.rs               # Peptide system management
│   ├── sampling.rs             # Sobol sequence generation
│   ├── forcefield.rs           # AMBER energy calculations
│   ├── minimizer.rs            # L-BFGS optimization
│   ├── parser.rs               # PDB I/O and parsing
│   └── main.rs                 # CLI interface
│
├── python/                     # Python bindings
│   └── src/
│       ├── lib.rs              # PyO3 interface
│       └── main.py             # OpenMM integration
│
├── library/                    # Force field parameters
├── examples/                   # Usage examples
├── docs/                       # Documentation
└── tests/                      # Test suite
```

---

## 📚 Documentation

| Resource | Description |
|----------|-------------|
| [**Tutorial**](docs/TUTORIAL.md) | Step-by-step guide from basics to advanced usage |
| [**API Reference**](docs/API_REFERENCE.md) | Complete function and class documentation |
| [**Python API**](docs/PYTHON_API.md) | Python-specific interface guide |
| [**System Module**](docs/system. md) | Peptide system management |
| [**Sampling Module**](docs/sampling.md) | Sobol sequence implementation |
| [**Force Fields**](docs/forcefield.md) | AMBER energy calculations |
| [**Minimizer**](docs/minimizer.md) | L-BFGS optimization details |

---

## 📊 Performance

DNCS is optimized for efficient peptide structure prediction:

| Peptide Size | Samples | Approximate Time* |
|--------------|---------|-------------------|
| 5 residues | 1,000 | ~2 seconds |
| 10 residues | 1,000 | ~5 seconds |
| 20 residues | 1,000 | ~15 seconds |
| 50 residues | 1,000 | ~60 seconds |

*Benchmarked on Intel i7-10700K, 8 cores, using `fold` method

---

## 🤝 Contributing

We welcome contributions! Please follow these guidelines:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5.  **Open** a Pull Request

### Reporting Issues

When reporting issues, please include:
- DNCS version (`dncs --version`)
- Operating system and architecture
- Minimal reproducible example
- Complete error messages and stack traces

---

## 📜 Citation

If you use DNCS in your research, please cite:

```bibtex
@article{dncs2024,
  title     = {Digital Nets Conformational Sampling for Peptide Structure Prediction},
  author    = {Rebairo, Abraham J. and Satheeshkumar, S. and Sam Paul, D.  and 
               Stephen, A. and {FiveAngstrom LLP}},
  journal   = {Physical Chemistry Chemical Physics},
  year      = {2024},
  doi       = {10. 1039/D4CP01891E},
  publisher = {Royal Society of Chemistry}
}
```

---

## 📄 License

This project is licensed under the **GNU General Public License v3. 0** - see the [LICENSE](LICENSE) file for details.

---

## 👥 Authors

**DNCS** is developed and maintained by:

- **Abraham Rebairo J.**
- **Satheeshkumar S.**
- **Sam Paul D.**
- **Stephen A.**
- **[FiveAngstrom LLP](https://fiveangstrom. com)**

---

<div align="center">

### 🧬 Happy Structure Prediction! ✨

*Bringing mathematical elegance to peptide modeling*

**[⬆ Back to Top](#-dncs)**

</div>
