# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DNCS (Digital Nets Conformational Sampling) is a molecular simulation library written in Rust with Python bindings. It provides enhanced conformational sampling using digital nets (Sobol sequences) for molecular dynamics simulations.

### Core Architecture

The codebase follows a modular architecture with these main components:

- **System** (`src/system.rs`): Core molecular system management, handles amino acid sequences, atom properties, and PDB I/O
- **Sampling** (`src/sampling.rs`): Implements Sobol sequence-based sampling with methods: Fold, Search, and Explore
- **Force Fields** (`src/forcefield.rs`): AMBER force field implementation supporting FF03, FF10, FF96, FF99SB, and FB15
- **Minimizer** (`src/minimizer.rs`): L-BFGS energy minimization algorithms
- **Parser** (`src/parser.rs`): PDB file parsing and molecular data structures
- **Python Bindings** (`python/src/lib.rs`): PyO3-based Python interface with OpenMM integration

### Key Design Patterns

- Uses `Arc<System>` for thread-safe shared molecular system access
- Parallel processing via `rayon` crate for sampling operations
- Configuration through TOML files (`dncs.toml` for Python interface)
- Force field parameters loaded from TOML/JSON files in `library/` directory

## Development Commands

### Rust Development
```bash
# Build the library
cargo build --release

# Run tests
cargo test

# Run specific test suite
./tests/run_tests.sh accuracy    # Forcefield accuracy tests
./tests/run_tests.sh components  # Force component tests
./tests/run_tests.sh quick       # Quick test subset

# Run examples
cargo run --example energy
cargo run --example sample
```

### Python Development
```bash
# Install dependencies and build Python bindings
just install

# Install only Python bindings (after deps installed)
just install-dncs

# Run Python interface
just run

# Manual installation
pip install -r python/requirements.txt
maturin develop --release -m python/Cargo.toml
```

### Build Tools
- **Just**: Primary build tool (`justfile`) for Python workflows
- **Maturin**: Used for building Python bindings
- **Cargo**: Standard Rust build system

## Testing Strategy

The project has comprehensive testing:

- `tests/forcefield_accuracy_tests.rs`: Energy calculation validation against reference values
- `tests/force_components_tests.rs`: Individual force component testing (bonds, angles, dihedrals, non-bonded)
- `tests/run_tests.sh`: Test runner script with multiple test categories

Use `./tests/run_tests.sh help` to see all available test options.

## Configuration Files

### For Rust Interface
- Use library files in `library/` directory (TOML format for force field parameters)
- Examples in `examples/` directory show typical usage patterns

### For Python Interface
Create `dncs.toml` with required fields:
- `moleculename`: Molecule identifier
- `sequence`: Amino acid sequence
- `forcefield`: Array of force field XML files
- `method`: Sampling method ("fold", "search", or "explore")
- `n_samples`: Number of samples to generate
- `temp`: Temperature for sampling

## File Structure Conventions

- Force field XML files in `library/ForceFields/`
- Configuration templates in root directory (`dncs.toml`)
- Results output to `Result/` directory by default
- Python integrator script at `python/src/integrator.py`

## Python-OpenMM Integration

The Python interface (`python/src/integrator.py`) integrates with OpenMM for:
- Langevin dynamics equilibration
- NPT/NVT molecular dynamics simulations
- Multi-step simulation workflows (sampling → minimization → equilibration → production)

When working with the Python interface, ensure OpenMM is properly installed and configured for the target device (CPU/CUDA/HIP).