use std::f64::consts::PI;



// Calculate the harmonic bond energy
//
fn harmonic_bond_energy(i: &Atom,j: &Atom) -> f64{
    let r = Self:distance(i,j);
    let r0 = r*0.5 //Equilibrium bond length
    let k_b= 1.00 //Force Constant
    0.5 * k_b * (r-r0).powi(2)
}
 
// Calculate the harmonic angle energy
fn harmonic_angle_energy(i: &Atom, j: &Atom, k: &Atom) -> f64{

    let theta = Self:calculate_bond_angle(i,j,k);
    let theta0 = theta*0.5; //equilibrium_angle
    let k_theta = 1.00 //Angular Force constant.
    let theta_rad = theta * PI / 180.0;   // Convert degrees to radians
    let theta0_rad = theta0 * PI / 180.0; // Convert degrees to radians
    0.5 * k_theta * (theta_rad - theta0_rad).powi(2)
}



// Calculate the periodic torsion energy
fn periodic_torsion_energy(i: &Atom, j: &Atom, k: &Atom, l: &Atom) -> f64{

    let phi = calculate_dihedral_angle(i,j,k,l);
    let phi_rad = phi * PI / 180.0;       // Convert degrees to radians
    let torsion = phi*0.5; //torsion phase and periodicity
    let gamma_rad = torsion.phase * PI / 180.0; // Convert degrees to radians

    0.5 * torsion.barrier_height * (1.0 + (torsion.periodicity as f64 * phi_rad - gamma_rad).cos())
}

fn periodic_torsion_energy(phi: f64, torsion: &Torsion) -> f64 {
    let phi_rad = phi * PI / 180.0;       // Convert degrees to radians
    let gamma_rad = torsion.phase * PI / 180.0; // Convert degrees to radians
    0.5 * torsion.barrier_height * (1.0 + (torsion.periodicity as f64 * phi_rad - gamma_rad).cos())
}



fn total_bonded_energy(atoms: &Vec<Atom>, bonds: &Vec<Bond>, angles: &Vec<Angle>, torsions: &Vec<Torsion>) -> f64 {
    let mut total_energy = 0.0;

    // HashMap for efficient atom lookup
    let atom_map: HashMap<usize, &Atom> = atoms.iter().map(|a| (a.id, a)).collect();

    // Calculate Bond Energies
    for bond in bonds {
        let a1 = match atom_map.get(&bond.atom1_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", bond.atom1_id);
                continue;
            }
        };

        let a2 = match atom_map.get(&bond.atom2_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", bond.atom2_id);
                continue;
            }
        };

        let r = calculate_distance(a1, a2);
        let energy = harmonic_bond_energy(r, bond.equilibrium_length, bond.force_constant);
        total_energy += energy;

        println!(
            "Bond {}-{}: r = {:.3} Å, r0 = {:.3} Å, k_b = {:.2} kcal/mol·Å², E = {:.3} kcal/mol",
            a1.name, a2.name, r, bond.equilibrium_length, bond.force_constant, energy
        );
    }

    // Calculate Angle Energies
    for angle in angles {
        let a1 = match atom_map.get(&angle.atom1_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", angle.atom1_id);
                continue;
            }
        };

        let a2 = match atom_map.get(&angle.atom2_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", angle.atom2_id);
                continue;
            }
        };

        let a3 = match atom_map.get(&angle.atom3_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", angle.atom3_id);
                continue;
            }
        };

        let theta = calculate_bond_angle(a1, a2, a3);
        let energy = harmonic_angle_energy(theta, angle.equilibrium_angle, angle.force_constant);
        total_energy += energy;

        println!(
            "Angle {}-{}-{}: theta = {:.3}°, theta0 = {:.3}°, k_theta = {:.2} kcal/mol·rad², E = {:.3} kcal/mol",
            a1.name, a2.name, a3.name, theta, angle.equilibrium_angle, angle.force_constant, energy
        );
    }

    // Calculate Torsion Energies
    for torsion in torsions {
        let a1 = match atom_map.get(&torsion.atom1_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", torsion.atom1_id);
                continue;
            }
        };

        let a2 = match atom_map.get(&torsion.atom2_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", torsion.atom2_id);
                continue;
            }
        };

        let a3 = match atom_map.get(&torsion.atom3_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", torsion.atom3_id);
                continue;
            }
        };

        let a4 = match atom_map.get(&torsion.atom4_id) {
            Some(atom) => atom,
            None => {
                eprintln!("Error: Atom ID {} not found.", torsion.atom4_id);
                continue;
            }
        };

        let phi = calculate_dihedral_angle(a1, a2, a3, a4);
        let energy = periodic_torsion_energy(phi, torsion);
        total_energy += energy;

        println!(
            "Torsion {}-{}-{}-{}: phi = {:.3}°, V_n = {:.2} kcal/mol, n = {}, gamma = {:.2}°, E = {:.3} kcal/mol",
            a1.name, a2.name, a3.name, a4.name, phi, torsion.barrier_height, torsion.periodicity, torsion.phase, energy
        );
    }

    println!("\nTotal Potential Energy: {:.3} kcal/mol", total_energy);

    total_energy
}

// Calculate Euclidean distance between two atoms
fn calculate_distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let dx = atom1.position[0] - atom2.position[0];
    let dy = atom1.position[1] - atom2.position[1];
    let dz = atom1.position[2] - atom2.position[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Calculate the bond angle (in degrees) between three atoms: A-B-C
fn calculate_bond_angle(atom_a: &Atom, atom_b: &Atom, atom_c: &Atom) -> f64 {
    // Vectors BA and BC
    let ba = [
        atom_a.position[0] - atom_b.position[0],
        atom_a.position[1] - atom_b.position[1],
        atom_a.position[2] - atom_b.position[2],
    ];
    let bc = [
        atom_c.position[0] - atom_b.position[0],
        atom_c.position[1] - atom_b.position[1],
        atom_c.position[2] - atom_b.position[2],
    ];
    
    // Dot product and magnitudes
    let dot_product = ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2];
    let mag_ba = (ba[0].powi(2) + ba[1].powi(2) + ba[2].powi(2)).sqrt();
    let mag_bc = (bc[0].powi(2) + bc[1].powi(2) + bc[2].powi(2)).sqrt();
    
    // Avoid division by zero
    if mag_ba == 0.0 || mag_bc == 0.0 {
        return 0.0;
    }
    
    // Calculate the angle in radians and then convert to degrees
    let cos_theta = dot_product / (mag_ba * mag_bc);
    let theta_rad = cos_theta.acos();
    let theta_deg = theta_rad * 180.0 / PI;
    
    theta_deg
}
/// Calculate the dihedral angle (in degrees) between four atoms: A-B-C-D
fn calculate_dihedral_angle(atom_a: &Atom, atom_b: &Atom, atom_c: &Atom, atom_d: &Atom) -> f64 {
    // Vectors BA, BC, and CD
    let ba = [
        atom_a.position[0] - atom_b.position[0],
        atom_a.position[1] - atom_b.position[1],
        atom_a.position[2] - atom_b.position[2],
    ];
    let bc = [
        atom_c.position[0] - atom_b.position[0],
        atom_c.position[1] - atom_b.position[1],
        atom_c.position[2] - atom_b.position[2],
    ];
    let cd = [
        atom_d.position[0] - atom_c.position[0],
        atom_d.position[1] - atom_c.position[1],
        atom_d.position[2] - atom_c.position[2],
    ];
    
    // Cross products
    let n1 = [
        ba[1]*bc[2] - ba[2]*bc[1],
        ba[2]*bc[0] - ba[0]*bc[2],
        ba[0]*bc[1] - ba[1]*bc[0],
    ];
    let n2 = [
        bc[1]*cd[2] - bc[2]*cd[1],
        bc[2]*cd[0] - bc[0]*cd[2],
        bc[0]*cd[1] - bc[1]*cd[0],
    ];
    
    // Normalize the vectors
    let mag_n1 = (n1[0].powi(2) + n1[1].powi(2) + n1[2].powi(2)).sqrt();
    let mag_n2 = (n2[0].powi(2) + n2[1].powi(2) + n2[2].powi(2)).sqrt();
    
    if mag_n1 == 0.0 || mag_n2 == 0.0 {
        return 0.0;
    }
    
    let n1_normalized = [n1[0]/mag_n1, n1[1]/mag_n1, n1[2]/mag_n1];
    let n2_normalized = [n2[0]/mag_n2, n2[1]/mag_n2, n2[2]/mag_n2];
    
    // Calculate the angle
    let cos_phi = n1_normalized[0]*n2_normalized[0] + n1_normalized[1]*n2_normalized[1] + n1_normalized[2]*n2_normalized[2];
    let sin_phi = (bc[0]*n1_normalized[1]*n2_normalized[2] +
                  bc[1]*n1_normalized[2]*n2_normalized[0] +
                  bc[2]*n1_normalized[0]*n2_normalized[1] -
                  bc[0]*n1_normalized[2]*n2_normalized[1] -
                  bc[1]*n1_normalized[0]*n2_normalized[2] -
                  bc[2]*n1_normalized[1]*n2_normalized[0]).abs();
    
    let phi_rad = cos_phi.acos();
    let phi_deg = if sin_phi >= 0.0 { phi_rad * 180.0 / PI } else { -phi_rad * 180.0 / PI };
    
    phi_deg
}

