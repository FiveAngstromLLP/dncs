# Forcefield Accuracy Test Suite

This directory contains comprehensive unit tests to validate the accuracy and reliability of the AMBER forcefield implementation in DNCS (Digital Nets Conformational Sampling).

## Test Files

### `forcefield_accuracy_tests.rs`
Main test suite for overall forcefield accuracy validation:

- **Energy Calculation Tests**: Validates consistency, scaling, and properties of energy calculations
- **Forcefield Variants**: Tests different AMBER forcefield versions (03, 96, 99SB, FB15)
- **System Size Scaling**: Verifies energy behavior with increasing peptide lengths
- **Charged Residue Tests**: Validates electrostatic interactions for charged amino acids
- **PDB Loading Tests**: Tests energy calculation from PDB structures
- **Performance Tests**: Benchmarks calculation speed and memory usage
- **Reproducibility Tests**: Ensures deterministic results across multiple runs
- **Integration Tests**: Complete workflow validation with various sequences

### `force_components_tests.rs`
Detailed validation of individual force components:

- **Harmonic Bond Forces**: Tests bond length and force constant parameters
- **Harmonic Angle Forces**: Validates angle parameters and calculations
- **Periodic Torsion Forces**: Tests dihedral angle parameters and periodicity
- **Non-bonded Interactions**: Validates Lennard-Jones and Coulomb parameters
- **Distance/Angle Calculations**: Tests geometric calculation accuracy
- **Parameter Bounds**: Validates physical reasonableness of all parameters
- **Numerical Precision**: Tests for numerical stability and precision
- **Parameter Coverage**: Ensures all amino acids have proper parameters

## Running the Tests

### Run All Tests
```bash
cargo test --test forcefield_accuracy_tests
cargo test --test force_components_tests
```

### Run Specific Test Categories
```bash
# Energy calculation tests
cargo test --test forcefield_accuracy_tests test_energy

# Force component tests
cargo test --test force_components_tests test_harmonic

# Performance tests
cargo test --test forcefield_accuracy_tests test_performance
```

### Run with Output
```bash
# See println! output from tests
cargo test --test forcefield_accuracy_tests -- --nocapture

# Run specific test with output
cargo test --test forcefield_accuracy_tests test_integration -- --nocapture
```

## Test Categories

### 1. Accuracy Tests
- Energy calculation consistency
- Cross-validation between forcefields
- Comparison with expected physical behavior
- Parameter validation against literature values

### 2. Reliability Tests
- Numerical stability and precision
- Reproducibility across runs
- Error handling for edge cases
- Memory safety and performance

### 3. Coverage Tests
- All amino acid types
- Different system sizes
- Various molecular conformations
- Multiple forcefield variants

### 4. Integration Tests
- Complete workflows from sequence to energy
- PDB file loading and processing
- System initialization and parameter assignment
- Multi-component energy calculations

## Expected Behavior

### Energy Ranges
- **Small peptides (1-3 residues)**: -1000 to +1000 kcal/mol
- **Medium peptides (4-10 residues)**: -5000 to +5000 kcal/mol
- **Large systems**: Depends on conformation and interactions

### Parameter Ranges
- **Bond lengths**: 0.5 - 3.0 Å
- **Bond force constants**: 0 - 10000 kcal/mol/Å²
- **Angles**: 60° - 180°
- **Angle force constants**: 0 - 1000 kcal/mol/rad²
- **Atomic masses**: 0.5 - 300 amu
- **LJ sigma**: 1.0 - 6.0 Å
- **LJ epsilon**: 0 - 1.0 kcal/mol
- **Partial charges**: -2.0 to +2.0 e

### Performance Expectations
- **Energy calculation time**: < 10 seconds for typical peptides
- **Memory usage**: Reasonable scaling with system size
- **Numerical precision**: Reproducible to machine precision

## Troubleshooting

### Common Issues

1. **Test Failures Due to Missing Files**
   - Ensure `library/ForceFields/test.pdb` exists for PDB tests
   - Check that all forcefield XML files are present

2. **Energy Values Outside Expected Ranges**
   - May indicate parameter assignment issues
   - Check atom typing and parameter lookup
   - Verify system initialization

3. **Numerical Precision Issues**
   - Could indicate floating-point instability
   - Check for division by zero or very small numbers
   - Verify unit conversions

4. **Performance Issues**
   - Large systems may take longer than expected
   - Check for inefficient algorithms or memory leaks
   - Consider parallel computation issues

### Debugging Tips

1. **Use verbose output**: Add `-- --nocapture` to see test output
2. **Run individual tests**: Isolate specific failing tests
3. **Check parameter assignment**: Verify atoms have proper types and parameters
4. **Validate input structures**: Ensure reasonable starting geometries

## Adding New Tests

When adding new tests:

1. **Follow naming conventions**: Use descriptive test function names
2. **Add documentation**: Explain what each test validates
3. **Use appropriate tolerances**: Set realistic comparison thresholds
4. **Test edge cases**: Include boundary conditions and error cases
5. **Validate assumptions**: Test underlying assumptions about the forcefield

## Test Data

### Reference Sequences
- `"AA"`: Simple dipeptide for basic tests
- `"AAA"`: Tripeptide for angle tests
- `"YGGFM"`: Mixed sequence with different residue types
- `"DEKR"`: Charged residues for electrostatic tests
- All 20 amino acids: Complete parameter coverage

### Test Structures
- `library/ForceFields/test.pdb`: Large protein structure for comprehensive testing
- Generated sequences: Systematic testing of all amino acid combinations

## Contributing

When contributing to the test suite:

1. Ensure tests are deterministic and reproducible
2. Use appropriate assertions and error messages
3. Test both success and failure cases
4. Document any new test categories or methodologies
5. Validate tests against known reference values when possible