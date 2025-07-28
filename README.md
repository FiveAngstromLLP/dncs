# DNCS Documentation

Welcome to the comprehensive documentation for DNCS (Digital Nets Conformational Sampling), a powerful molecular simulation library for enhanced conformational sampling using digital nets.

## Quick Navigation

### 📚 **Getting Started**
- [**Tutorial**](TUTORIAL.md) - Complete guide from basics to advanced usage
- [**API Reference**](API_REFERENCE.md) - Comprehensive API documentation
- [**Python API**](PYTHON_API.md) - Detailed Python interface documentation

### 🔧 **Core Components**
- [**System**](system.md) - Molecular system management
- [**Sampling**](sampling.md) - Conformational sampling methods
- [**Force Fields**](forcefield.md) - Energy calculations and force fields
- [**Minimizer**](minimizer.md) - Energy minimization
- [**Parser**](parser.md) - File parsing and data structures
- [**CLI**](cli.md) - Command-line interface

## Overview

DNCS is designed for molecular simulation researchers who need:

- **Enhanced Sampling**: Uniform exploration of conformational space using Sobol sequences
- **Multiple Interfaces**: Both Rust and Python APIs for flexibility
- **Force Field Support**: Multiple AMBER force fields (03, 10, 96, 99SB, FB15)
- **Integration**: Seamless integration with OpenMM for molecular dynamics
- **Performance**: Parallel processing and optimized algorithms

## Documentation Structure

### 🎯 **For New Users**

1. **Start Here**: [Tutorial](TUTORIAL.md)
   - Installation and setup
   - Basic concepts and theory
   - Step-by-step examples
   - Best practices

2. **Quick Reference**: [API Reference](API_REFERENCE.md)
   - Complete API documentation
   - Code examples for all functions
   - Configuration options
   - Error handling

### 🔬 **For Python Users**

- [**Python API Documentation**](PYTHON_API.md)
  - Complete Python interface
  - Integration examples
  - Performance tips
  - Advanced usage patterns

### ⚙️ **For Advanced Users**

- **Component Documentation**:
  - [System](system.md) - Core molecular system class
  - [Sampling](sampling.md) - Sobol sequences and sampling methods
  - [Force Fields](forcefield.md) - AMBER force field implementation
  - [Minimizer](minimizer.md) - L-BFGS energy minimization
  - [Parser](parser.md) - PDB parsing and data structures

## Key Features

### 🧬 **Molecular Systems**
- Support for amino acid sequences
- PDB file import/export
- Multiple force field support
- Comprehensive atom property management

### 📊 **Conformational Sampling**
- **Sobol Sequences**: Low-discrepancy sampling for uniform coverage
- **Multiple Methods**: Fold, Search, and Explore sampling strategies
- **Temperature Control**: Boltzmann weighting at different temperatures
- **Parallel Processing**: Efficient multi-core utilization

### ⚡ **Performance**
- **Rust Core**: High-performance implementation
- **Python Bindings**: Easy-to-use interface
- **Memory Efficient**: Optimized memory management
- **Scalable**: Handles large molecular systems

### 🔗 **Integration**
- **OpenMM**: Molecular dynamics simulation
- **Standard Formats**: PDB input/output
- **Flexible Configuration**: TOML configuration files

## Quick Start Examples

### Rust
```rust
use std::sync::Arc;
use libdncs::*;

fn main() {
    // Create molecular system
    let mut system = System::new("YGGFM", FF::AmberFB15.init());
    system.init_parameters();
    
    // Calculate energy
    let amber = Amber::new(Arc::new(system.clone()));
    let energy = amber.energy();
    println!("Energy: {:.2} kcal/mol", energy);
    
    // Perform sampling
    let mut sampler = Sampler::new(
        Arc::new(system),
        Method::Fold,
        "results".to_string()
    );
    sampler.sample(1000, 300.0);
}
```

### Python
```python
import dncs

# Create polymer system
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")

# Calculate energy
energy = polymer.getEnergy()
print(f"Energy: {energy:.2f} kcal/mol")

# Generate samples
sampler = dncs.SobolSampler(
    polymer, 1000, "fold", 300.0, "results"
)

# Save structure
polymer.toPDB("structure.pdb")
```

## Installation

### Prerequisites
- Rust toolchain (1.70.0+)
- Python 3.7+ (for Python bindings)
- Maturin (for building Python bindings)

### Quick Install
```bash
# Clone repository
git clone <repository-url>
cd dncs

# Install using just
just install

# Or install manually
pip install -r requirements.txt
pip install maturin
maturin develop --release
```

## Documentation Conventions

### Code Examples
- All examples are tested and working
- Rust examples use `libdncs::*`
- Python examples use `import dncs`
- Error handling is included where appropriate

### API Documentation
- **Parameters**: Type and description for all parameters
- **Returns**: Return type and description
- **Examples**: Working code examples for each function
- **Notes**: Important usage notes and limitations

### File Organization
```
docs/
├── README.md           # This file - documentation index
├── TUTORIAL.md         # Complete tutorial guide
├── API_REFERENCE.md    # Comprehensive API reference
├── PYTHON_API.md       # Python-specific documentation
├── system.md           # Core system documentation
├── sampling.md         # Sampling methods
├── forcefield.md       # Force field implementation
├── minimizer.md        # Energy minimization
├── parser.md           # File parsing
└── cli.md             # Command-line interface
```

## Support and Contributing

### Getting Help
1. Check the [Tutorial](TUTORIAL.md) for step-by-step guidance
2. Consult the [API Reference](API_REFERENCE.md) for detailed function documentation
3. Review examples in component documentation
4. Check the troubleshooting section in the tutorial

### Reporting Issues
When reporting issues, please include:
- DNCS version
- Operating system
- Minimal reproducible example
- Error messages and stack traces
- Expected vs. actual behavior

### Contributing
- Follow existing documentation style
- Include working examples
- Update relevant documentation sections
- Test all code examples

## Version Information

- **Current Version**: 1.1.2
- **Minimum Rust**: 1.70.0
- **Python Support**: 3.7+
- **License**: GPL v3

## Citation

If you use DNCS in your research, please cite:

```
Digital Nets Conformational Sampling (DNCS)
DOI: 10.1039/D4CP01891E
```

---

**Happy Sampling!** 🧬✨

For questions or support, please refer to the documentation sections above or check the repository issues.