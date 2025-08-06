use libdncs::*;

fn main() {
    println!("Detailed energy analysis of clash structure (Sample 854)...");

    let mut sys = System::from_pdb("Result/6RRO/sample_0854.pdb", FF::AmberFB15.init());
    sys.init_parameters();

    // Print first few atoms to see charges and positions
    println!("\nFirst 10 atoms:");
    for i in 0..10.min(sys.particles.len()) {
        let atom = &sys.particles[i];
        println!(
            "Atom {}: {} charge={:.3} pos=({:.3}, {:.3}, {:.3})",
            atom.serial,
            atom.name,
            atom.charge,
            atom.position[0],
            atom.position[1],
            atom.position[2]
        );
    }

    // Look for very short distances
    println!("\nSearching for very short distances (<0.15 nm):");
    let mut short_distances = Vec::new();

    for i in 0..sys.particles.len() {
        for j in (i + 1)..sys.particles.len() {
            let atom_i = &sys.particles[i];
            let atom_j = &sys.particles[j];

            let dx = atom_j.position[0] - atom_i.position[0];
            let dy = atom_j.position[1] - atom_i.position[1];
            let dz = atom_j.position[2] - atom_i.position[2];
            let distance = (dx * dx + dy * dy + dz * dz).sqrt() * 0.1; // Convert to nm

            if distance < 0.15 {
                // Calculate individual energies for this pair
                let lj_energy = if distance > 0.001 {
                    let sigma = (atom_i.sigma + atom_j.sigma) / 2.0;
                    let epsilon = (atom_i.epsilon * atom_j.epsilon).sqrt();
                    4.0 * epsilon * ((sigma / distance).powi(12) - (sigma / distance).powi(6))
                } else {
                    f64::INFINITY
                };

                let elec_energy = if distance > 0.001 {
                    138.93545727242866 * atom_i.charge * atom_j.charge / distance
                } else {
                    if atom_i.charge * atom_j.charge > 0.0 {
                        f64::INFINITY
                    } else {
                        f64::NEG_INFINITY
                    }
                };

                short_distances.push((i, j, distance, lj_energy, elec_energy));
            }
        }
    }

    // Sort by distance
    short_distances.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

    println!(
        "Found {} atom pairs with distance < 0.15 nm:",
        short_distances.len()
    );
    for (i, (idx1, idx2, dist, lj, elec)) in short_distances.iter().take(10).enumerate() {
        let atom1 = &sys.particles[*idx1];
        let atom2 = &sys.particles[*idx2];
        println!(
            "{}. Atoms {}({}) - {}({}) dist={:.4}nm LJ={:.1} Elec={:.1} q1={:.3} q2={:.3}",
            i + 1,
            atom1.serial,
            atom1.name,
            atom2.serial,
            atom2.name,
            dist,
            lj,
            elec,
            atom1.charge,
            atom2.charge
        );
    }

    let total_lj: f64 = short_distances
        .iter()
        .map(|(_, _, _, lj, _)| if lj.is_finite() { *lj } else { 1000.0 })
        .sum();
    let total_elec: f64 = short_distances
        .iter()
        .map(|(_, _, _, _, elec)| {
            if elec.is_finite() {
                *elec
            } else {
                if *elec > 0.0 {
                    1000.0
                } else {
                    -1000.0
                }
            }
        })
        .sum();

    println!("\nFrom short distances only:");
    println!("Total LJ contribution: {:.1} kJ/mol", total_lj);
    println!("Total Electrostatic contribution: {:.1} kJ/mol", total_elec);
    println!("Net from clashes: {:.1} kJ/mol", total_lj + total_elec);
}
