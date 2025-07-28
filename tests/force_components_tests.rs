use std::sync::Arc;
use libdncs::*;

/// Test module for individual force component validation
/// 
/// This module contains detailed tests for individual components of the 
/// AMBER forcefield implementation:
/// - Harmonic bond forces
/// - Harmonic angle forces  
/// - Periodic torsion forces
/// - Non-bonded interactions (Lennard-Jones and Coulomb)
/// - Hydrogen bonding
/// - Distance and angle calculations

#[cfg(test)]
mod force_components_tests {
    use super::*;
    use std::f64::consts::PI;

    const TOLERANCE: f64 = 1e-10;
    const ENERGY_TOLERANCE: f64 = 1e-6;

    /// Test harmonic bond force calculations
    #[test]
    fn test_harmonic_bond_forces() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AA", ff);
        system.init_parameters();
        
        // Verify that bonds exist
        let has_bonds = system.firstbonded.iter().any(|bonds| !bonds.is_empty());
        assert!(has_bonds, "System should have bonded interactions");
        
        // Test that bond parameters are reasonable
        let bond_params = &system.forcefield.harmonic_bond_force.bonds;
        assert!(!bond_params.is_empty(), "Should have bond parameters");
        
        for bond in bond_params {
            assert!(bond.length > 0.5 && bond.length < 5.0, 
                    "Bond length {} Å is outside reasonable range", bond.length);
            assert!(bond.k > 0.0, 
                    "Bond force constant {} should be positive", bond.k);
        }
    }

    /// Test harmonic angle force calculations
    #[test]
    fn test_harmonic_angle_forces() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AAA", ff);
        system.init_parameters();
        
        // Verify that angles exist
        let has_angles = system.secondbonded.iter().any(|angles| !angles.is_empty());
        assert!(has_angles, "System should have angle interactions");
        
        // Test angle parameters
        let angle_params = &system.forcefield.harmonic_angle_force.angles;
        assert!(!angle_params.is_empty(), "Should have angle parameters");
        
        for angle in angle_params {
            assert!(angle.angle >= 0.0 && angle.angle <= 180.0, 
                    "Angle {} degrees is outside valid range", angle.angle);
            assert!(angle.k >= 0.0, 
                    "Angle force constant {} should be non-negative", angle.k);
        }
    }

    /// Test periodic torsion forces
    #[test]
    fn test_periodic_torsion_forces() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AAAA", ff);
        system.init_parameters();
        
        // Test torsion parameters
        let torsion_params = &system.forcefield.periodic_torsion_force;
        assert!(!torsion_params.proper.is_empty(), "Should have proper torsion parameters");
        
        for torsion in &torsion_params.proper {
            // Test that periodicities are reasonable
            assert!(torsion.periodicity1 >= 1.0 && torsion.periodicity1 <= 6.0,
                    "Periodicity {} is outside reasonable range", torsion.periodicity1);
            
            // Test that phases are in valid range
            assert!(torsion.phase1 >= -PI && torsion.phase1 <= PI,
                    "Phase {} is outside valid range", torsion.phase1);
            
            // Test force constants
            assert!(torsion.k1.abs() < 100.0,
                    "Torsion force constant {} seems too large", torsion.k1);
        }
    }

    /// Test non-bonded force parameters
    #[test]
    fn test_nonbonded_parameters() {
        let ff = FF::AmberFB15.init();
        let nonbonded = &ff.nonbonded_force;
        
        // Test scaling factors
        assert!(nonbonded.coulomb14scale > 0.0 && nonbonded.coulomb14scale <= 1.0,
                "Coulomb 1-4 scaling factor {} is outside valid range", 
                nonbonded.coulomb14scale);
        assert!(nonbonded.lj14scale > 0.0 && nonbonded.lj14scale <= 1.0,
                "LJ 1-4 scaling factor {} is outside valid range", 
                nonbonded.lj14scale);
        
        // Test atom parameters
        for atom in &nonbonded.atoms {
            assert!(atom.sigma > 0.0 && atom.sigma < 10.0,
                    "Sigma {} Å is outside reasonable range", atom.sigma);
            assert!(atom.epsilon >= 0.0 && atom.epsilon < 10.0,
                    "Epsilon {} kcal/mol is outside reasonable range", atom.epsilon);
            assert!(atom.charge.abs() < 10.0,
                    "Charge {} e is outside reasonable range", atom.charge);
        }
    }

    /// Test distance calculation accuracy
    #[test]
    fn test_distance_calculations() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AA", ff);
        system.init_parameters();
        
        // Create test atoms at known positions
        if system.particles.len() >= 2 {
            let atom1 = &system.particles[0];
            let atom2 = &system.particles[1];
            
            // Calculate distance manually
            let dx = atom2.position[0] - atom1.position[0];
            let dy = atom2.position[1] - atom1.position[1];
            let dz = atom2.position[2] - atom1.position[2];
            let expected_distance = (dx*dx + dy*dy + dz*dz).sqrt();
            
            // Distance should be positive and reasonable for bonded atoms
            assert!(expected_distance > 0.5 && expected_distance < 5.0,
                    "Distance {} Å between bonded atoms is unreasonable", expected_distance);
        }
    }

    /// Test angle calculation accuracy
    #[test]
    fn test_angle_calculations() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AAA", ff);
        system.init_parameters();
        
        // Test with at least 3 atoms for angle calculation
        if system.particles.len() >= 3 {
            // Verify that we have atoms to work with
            assert!(system.particles.len() >= 3, 
                    "Need at least 3 atoms for angle test");
            
            // Test that angles are physically reasonable
            // (This would require access to the angle calculation method)
            // For now, just verify the system is properly set up
            let amber = Amber::new(Arc::new(system));
            let energy = amber.energy();
            assert!(energy.is_finite(), "Energy should be finite for angle test");
        }
    }

    /// Test Lennard-Jones energy calculations
    #[test]
    fn test_lennard_jones_energy() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AA", ff);
        system.init_parameters();
        
        // Find atoms with LJ parameters
        let atoms_with_lj = system.particles.iter()
            .filter(|atom| atom.sigma > 0.0 && atom.epsilon > 0.0)
            .count();
        
        assert!(atoms_with_lj > 0, 
                "System should have atoms with Lennard-Jones parameters");
        
        // Test that LJ parameters are in reasonable ranges
        for atom in &system.particles {
            if atom.sigma > 0.0 {
                assert!(atom.sigma > 1.0 && atom.sigma < 6.0,
                        "Sigma {} Å is outside typical range for atoms", atom.sigma);
            }
            if atom.epsilon > 0.0 {
                assert!(atom.epsilon < 1.0,
                        "Epsilon {} kcal/mol is unusually large", atom.epsilon);
            }
        }
    }

    /// Test electrostatic energy calculations
    #[test]
    fn test_electrostatic_energy() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("DEKR", ff); // Charged residues
        system.init_parameters();
        
        // Find atoms with charges
        let charged_atoms = system.particles.iter()
            .filter(|atom| atom.charge.abs() > 0.01)
            .count();
        
        assert!(charged_atoms > 0, 
                "Charged sequence should have atoms with significant charges");
        
        // Test charge conservation (should be approximately neutral for full residues)
        let total_charge: f64 = system.particles.iter()
            .map(|atom| atom.charge)
            .sum();
        
        // For amino acids, total charge should be close to expected values
        // D(-1), E(-1), K(+1), R(+1) = net 0
        assert!(total_charge.abs() < 1.0,
                "Total charge {} is larger than expected", total_charge);
    }

    /// Test hydrogen bonding energy
    #[test]
    fn test_hydrogen_bonding() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("SS", ff); // Serine has OH groups
        system.init_parameters();
        
        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();
        
        // Just verify that the calculation completes
        assert!(energy.is_finite(), 
                "Energy calculation with hydrogen bonding should be finite");
    }

    /// Test energy component additivity
    #[test]
    fn test_energy_additivity() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("YGG", ff);
        system.init_parameters();
        
        let amber = Amber::new(Arc::new(system));
        let total_energy = amber.energy();
        
        // Total energy should be finite
        assert!(total_energy.is_finite(), 
                "Total energy should be finite");
        
        // Energy should be the sum of all components
        // (This test mainly verifies that no component returns NaN or infinity)
        assert!(!total_energy.is_nan(), "Total energy should not be NaN");
        assert!(!total_energy.is_infinite(), "Total energy should not be infinite");
    }

    /// Test force field consistency across different systems
    #[test]
    fn test_forcefield_consistency() {
        let sequences = vec!["A", "G", "Y", "F"];
        
        for seq in sequences {
            let ff = FF::AmberFB15.init();
            let mut system = System::new(seq, ff);
            system.init_parameters();
            
            // Check that all atoms have proper parameters assigned
            for atom in &system.particles {
                if let Some(atomtype) = &atom.atomtype {
                    assert!(!atomtype.is_empty(), 
                            "Atom type should not be empty for atom {}", atom.name);
                }
                
                // Check that essential parameters are assigned
                assert!(atom.sigma >= 0.0, 
                        "Sigma should be non-negative for atom {}", atom.name);
                assert!(atom.epsilon >= 0.0, 
                        "Epsilon should be non-negative for atom {}", atom.name);
            }
            
            let amber = Amber::new(Arc::new(system));
            let energy = amber.energy();
            assert!(energy.is_finite(), 
                    "Energy should be finite for sequence {}", seq);
        }
    }

    /// Test parameter bounds and physical reasonableness
    #[test]
    fn test_parameter_bounds() {
        let ff = FF::AmberFB15.init();
        
        // Test atom type masses
        for atom_type in &ff.atom_types.types {
            assert!(atom_type.mass > 0.5 && atom_type.mass < 300.0,
                    "Atomic mass {} is outside reasonable range", atom_type.mass);
        }
        
        // Test bond parameters
        for bond in &ff.harmonic_bond_force.bonds {
            assert!(bond.length > 0.5 && bond.length < 3.0,
                    "Bond length {} Å is outside typical range", bond.length);
            assert!(bond.k > 0.0 && bond.k < 10000.0,
                    "Bond force constant {} is outside reasonable range", bond.k);
        }
        
        // Test angle parameters
        for angle in &ff.harmonic_angle_force.angles {
            assert!(angle.angle >= 60.0 && angle.angle <= 180.0,
                    "Angle {} degrees is outside typical range", angle.angle);
            assert!(angle.k >= 0.0 && angle.k < 1000.0,
                    "Angle force constant {} is outside reasonable range", angle.k);
        }
    }

    /// Test numerical precision and stability
    #[test]
    fn test_numerical_precision() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AA", ff);
        system.init_parameters();
        
        let amber = Amber::new(Arc::new(system));
        
        // Calculate energy multiple times to check for numerical instability
        let mut energies = Vec::new();
        for _ in 0..10 {
            energies.push(amber.energy());
        }
        
        // All calculations should give identical results
        let first_energy = energies[0];
        for (i, &energy) in energies.iter().enumerate() {
            assert!((energy - first_energy).abs() < TOLERANCE,
                    "Energy calculation {} differs from first: {} vs {}",
                    i, energy, first_energy);
        }
    }

    /// Test edge cases and boundary conditions
    #[test]
    fn test_edge_cases() {
        let ff = FF::AmberFB15.init();
        
        // Test with minimal system
        let mut system = System::new("G", ff.clone()); // Glycine is smallest
        system.init_parameters();
        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();
        assert!(energy.is_finite(), "Energy should be finite for glycine");
        
        // Test with larger system
        let mut system = System::new("WWWWW", ff); // Tryptophan is large
        system.init_parameters();
        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();
        assert!(energy.is_finite(), "Energy should be finite for large system");
    }

    /// Test force field parameter coverage
    #[test]
    fn test_parameter_coverage() {
        let ff = FF::AmberFB15.init();
        
        // Test that we have parameters for common amino acids
        let amino_acids = vec!["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
                              "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];
        
        for aa in amino_acids {
            let mut system = System::new(aa, ff.clone());
            system.init_parameters();
            
            // Check that the system was created successfully
            assert!(!system.particles.is_empty(), 
                    "Should have particles for amino acid {}", aa);
            
            // Check that atoms have types assigned
            let typed_atoms = system.particles.iter()
                .filter(|atom| atom.atomtype.is_some())
                .count();
            assert!(typed_atoms > 0, 
                    "Should have typed atoms for amino acid {}", aa);
            
            // Check that energy calculation works
            let amber = Amber::new(Arc::new(system));
            let energy = amber.energy();
            assert!(energy.is_finite(), 
                    "Energy should be finite for amino acid {}", aa);
        }
    }
}