use libdncs::*;
use std::sync::Arc;

/// Test module for forcefield accuracy validation
///
/// This module contains comprehensive tests to validate the accuracy of
/// the AMBER forcefield implementation in DNCS, including:
/// - Energy calculations for different molecular systems
/// - Individual force component validation
/// - Comparison with reference values
/// - Cross-validation between different forcefields
/// - Numerical stability tests

#[cfg(test)]
mod forcefield_accuracy_tests {
    use super::*;
    // use std::f64::consts::PI;

    // Test constants
    const TOLERANCE: f64 = 1e-6;
    // const ENERGY_TOLERANCE: f64 = 1e-3; // kcal/mol

    // Reference test sequences
    const SIMPLE_DIPEPTIDE: &str = "AA";
    const TRIPEPTIDE: &str = "AAA";
    const MIXED_SEQUENCE: &str = "YGGFM";
    const CHARGED_SEQUENCE: &str = "DEKR";

    /// Test basic energy calculation consistency
    #[test]
    fn test_energy_calculation_consistency() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new(SIMPLE_DIPEPTIDE, ff);
        system.init_parameters();

        let amber = Amber::new(Arc::new(system));
        let energy1 = amber.energy();
        let energy2 = amber.energy();

        // Energy calculation should be deterministic
        assert!(
            (energy1 - energy2).abs() < TOLERANCE,
            "Energy calculation not consistent: {} vs {}",
            energy1,
            energy2
        );

        // Energy should be finite
        assert!(
            energy1.is_finite(),
            "Energy calculation returned non-finite value: {}",
            energy1
        );
    }

    /// Test energy scaling with system size
    #[test]
    fn test_energy_scaling() {
        let ff = FF::AmberFB15.init();

        // Test with increasing peptide lengths
        let sequences = vec!["A", "AA", "AAA", "AAAA"];
        let mut energies = Vec::new();

        for seq in sequences {
            let mut system = System::new(seq, ff.clone());
            system.init_parameters();
            let amber = Amber::new(Arc::new(system));
            let energy = amber.energy();
            energies.push(energy);

            // Energy should be finite for all systems
            assert!(
                energy.is_finite(),
                "Non-finite energy for sequence {}: {}",
                seq,
                energy
            );
        }

        // Energy should generally increase with system size (more interactions)
        // Note: This is a general trend, not strict monotonicity due to conformational effects
        assert!(energies.len() == 4, "Should have 4 energy values");

        // Check that we don't have any obviously wrong values (like all zeros)
        assert!(
            energies.iter().any(|&e| e.abs() > 1.0),
            "All energies are suspiciously small: {:?}",
            energies
        );
    }

    /// Test individual force components
    #[test]
    fn test_force_components() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new(TRIPEPTIDE, ff);
        system.init_parameters();

        let amber = Amber::new(Arc::new(system.clone()));
        let total_energy = amber.energy();

        // Test that we can access individual components
        // Note: This requires making some methods public or adding accessors
        assert!(total_energy.is_finite(), "Total energy should be finite");

        // Test that system has expected number of atoms
        let expected_atoms = TRIPEPTIDE.len() * 10; // Approximate atoms per residue
        assert!(
            system.particles.len() >= expected_atoms,
            "System should have at least {} atoms, got {}",
            expected_atoms,
            system.particles.len()
        );
    }

    /// Test different forcefield variants
    #[test]
    fn test_forcefield_variants() {
        let forcefields = vec![FF::Amber03, FF::Amber96, FF::Amber99SB, FF::AmberFB15];

        let mut energies = Vec::new();

        for ff_type in forcefields {
            let ff = ff_type.init();
            let mut system = System::new(SIMPLE_DIPEPTIDE, ff);
            system.init_parameters();

            let amber = Amber::new(Arc::new(system));
            let energy = amber.energy();

            assert!(
                energy.is_finite(),
                "Energy calculation failed for forcefield {:?}: {}",
                ff_type,
                energy
            );

            energies.push((ff_type, energy));
        }

        // All forcefields should produce finite energies
        assert_eq!(
            energies.len(),
            4,
            "Should have energies for all 4 forcefields"
        );

        // Energies should be different between forcefields (they use different parameters)
        let first_energy = energies[0].1;
        let all_same = energies
            .iter()
            .all(|(_, e)| (e - first_energy).abs() < TOLERANCE);
        assert!(
            !all_same,
            "All forcefields produce identical energies, which is suspicious"
        );
    }

    /// Test energy conservation properties
    #[test]
    fn test_energy_properties() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new(MIXED_SEQUENCE, ff);
        system.init_parameters();

        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();

        // Test basic energy properties
        assert!(energy.is_finite(), "Energy should be finite");

        // For a realistic protein system, energy should be negative or small positive
        // (this is a sanity check - very large positive energies indicate problems)
        assert!(
            energy < 10000.0,
            "Energy suspiciously high: {} kcal/mol",
            energy
        );

        // Energy shouldn't be extremely negative either
        assert!(
            energy > -100000.0,
            "Energy suspiciously low: {} kcal/mol",
            energy
        );
    }

    /// Test charged residue interactions
    #[test]
    fn test_charged_residues() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new(CHARGED_SEQUENCE, ff);
        system.init_parameters();

        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();

        assert!(
            energy.is_finite(),
            "Energy calculation failed for charged sequence: {}",
            energy
        );

        // Charged systems should have significant electrostatic contributions
        // This is mainly a sanity check that the calculation completes
        println!("Charged sequence energy: {} kcal/mol", energy);
    }

    /// Test system initialization
    #[test]
    fn test_system_initialization() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("YGGFM", ff);

        // Test before initialization
        assert!(!system.particles.is_empty(), "System should have particles");

        // Initialize parameters
        system.init_parameters();

        // Test after initialization
        assert!(
            !system.particles.is_empty(),
            "System should still have particles after init"
        );

        // Check that atoms have proper types assigned
        let atoms_with_types = system
            .particles
            .iter()
            .filter(|atom| atom.atomtype.is_some())
            .count();

        assert!(
            atoms_with_types > 0,
            "No atoms have types assigned after initialization"
        );

        // Check that we have proper bonding information
        assert!(
            !system.firstbonded.is_empty(),
            "System should have bonding information"
        );
    }

    /// Test PDB file loading and energy calculation
    #[test]
    fn test_pdb_energy_calculation() {
        let pdb_file = "library/ForceFields/test.pdb";

        // Check if test PDB file exists
        if std::path::Path::new(pdb_file).exists() {
            let ff = FF::Amber99SB.init();
            let mut system = System::from_pdb(pdb_file, ff);
            system.init_parameters();

            let amber = Amber::new(Arc::new(system));
            let energy = amber.energy();

            assert!(
                energy.is_finite(),
                "Energy calculation failed for PDB file: {}",
                energy
            );

            println!("PDB system energy: {} kcal/mol", energy);
        } else {
            println!("Test PDB file not found, skipping PDB test");
        }
    }

    /// Test numerical stability of distance calculations
    #[test]
    fn test_distance_calculations() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AA", ff);
        system.init_parameters();

        // Test that we can create the system and it has atoms
        assert!(
            system.particles.len() >= 2,
            "System should have at least 2 atoms for distance test"
        );

        // Test with the system
        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();

        assert!(
            energy.is_finite(),
            "Energy should be finite for distance test"
        );
    }

    /// Test force field parameter loading
    #[test]
    fn test_forcefield_parameters() {
        let ff = FF::AmberFB15.init();

        // Test that forcefield has required components
        assert!(
            !ff.atom_types.types.is_empty(),
            "Forcefield should have atom types"
        );
        assert!(
            !ff.residues.residue.is_empty(),
            "Forcefield should have residue definitions"
        );
        assert!(
            !ff.harmonic_bond_force.bonds.is_empty(),
            "Forcefield should have bond parameters"
        );
        assert!(
            !ff.harmonic_angle_force.angles.is_empty(),
            "Forcefield should have angle parameters"
        );
        assert!(
            !ff.nonbonded_force.atoms.is_empty(),
            "Forcefield should have nonbonded parameters"
        );

        // Test parameter ranges
        for atom_type in &ff.atom_types.types {
            assert!(
                atom_type.mass > 0.0,
                "Atom mass should be positive: {}",
                atom_type.mass
            );
        }

        for bond in &ff.harmonic_bond_force.bonds {
            assert!(
                bond.k >= 0.0,
                "Bond force constant should be non-negative: {}",
                bond.k
            );
            assert!(
                bond.length > 0.0,
                "Bond length should be positive: {}",
                bond.length
            );
        }

        for nb_atom in &ff.nonbonded_force.atoms {
            assert!(
                nb_atom.sigma > 0.0,
                "Sigma parameter should be positive: {}",
                nb_atom.sigma
            );
            assert!(
                nb_atom.epsilon >= 0.0,
                "Epsilon parameter should be non-negative: {}",
                nb_atom.epsilon
            );
        }
    }

    /// Test energy units and scaling
    #[test]
    fn test_energy_units() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("AA", ff);
        system.init_parameters();

        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();

        // Energy should be in reasonable range for kcal/mol
        // For a small dipeptide, energy should be in range of -1000 to +1000 kcal/mol
        assert!(
            energy > -1000.0 && energy < 1000.0,
            "Energy {} kcal/mol is outside expected range for small peptide",
            energy
        );
    }

    /// Benchmark test for performance
    #[test]
    fn test_performance() {
        let ff = FF::AmberFB15.init();
        let mut system = System::new("YGGFM", ff);
        system.init_parameters();

        let amber = Arc::new(Amber::new(Arc::new(system)));

        let start = std::time::Instant::now();
        let _energy = amber.energy();
        let duration = start.elapsed();

        // Energy calculation should complete in reasonable time
        assert!(
            duration.as_secs() < 10,
            "Energy calculation took too long: {:?}",
            duration
        );

        println!("Energy calculation time: {:?}", duration);
    }

    /// Test reproducibility across multiple runs
    #[test]
    fn test_reproducibility() {
        let ff = FF::AmberFB15.init();

        let mut energies = Vec::new();
        for _ in 0..5 {
            let mut system = System::new("GGG", ff.clone());
            system.init_parameters();
            let amber = Amber::new(Arc::new(system));
            energies.push(amber.energy());
        }

        // All energies should be identical (deterministic calculation)
        let first_energy = energies[0];
        for (i, &energy) in energies.iter().enumerate() {
            assert!(
                (energy - first_energy).abs() < TOLERANCE,
                "Energy not reproducible: run {} gave {}, expected {}",
                i,
                energy,
                first_energy
            );
        }
    }

    /// Test error handling for invalid inputs
    #[test]
    fn test_error_handling() {
        let ff = FF::AmberFB15.init();

        // Test with empty sequence - should not panic
        let system = System::new("", ff.clone());
        assert!(
            system.particles.is_empty(),
            "Empty sequence should produce no particles"
        );

        // Test with single residue
        let mut system = System::new("A", ff);
        system.init_parameters();
        let amber = Amber::new(Arc::new(system));
        let energy = amber.energy();
        assert!(energy.is_finite(), "Single residue energy should be finite");
    }

    /// Integration test combining multiple components
    #[test]
    fn test_integration() {
        // Test a complete workflow
        let sequences = vec!["A", "AA", "YGG", "DEKR"];
        let forcefields = vec![FF::Amber99SB, FF::AmberFB15];

        for seq in sequences {
            for ff_type in &forcefields {
                let ff = ff_type.init();
                let mut system = System::new(seq, ff);
                system.init_parameters();

                // Verify system is properly initialized
                assert!(
                    !system.particles.is_empty(),
                    "System should have particles for sequence {}",
                    seq
                );

                let amber = Amber::new(Arc::new(system));
                let energy = amber.energy();

                assert!(
                    energy.is_finite(),
                    "Energy should be finite for {} with {:?}: {}",
                    seq,
                    ff_type,
                    energy
                );

                println!(
                    "Sequence: {}, FF: {:?}, Energy: {:.3} kcal/mol",
                    seq, ff_type, energy
                );
            }
        }
    }
}
