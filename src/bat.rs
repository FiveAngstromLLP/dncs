// src/main.rs

use std::f64::consts::PI;


#[derive(Debug, Clone)]
struct Atom {
    id: usize,           // Unique identifier for the atom
    name: String,        // Atom name/type
    position: [f64; 3],  // x, y, z coordinates in Angstroms
}

#[derive(Debug, Clone)]
struct Bond {
    atom1_id: usize,                   // ID of the first atom
    atom2_id: usize,                   // ID of the second atom
    equilibrium_length: f64,           // Equilibrium bond length (r0) in Å
    force_constant: f64,               // Force constant (k_b) in kcal/mol·Å²
}
#[derive(Debug, Clone)]
struct Angle {
    atom1_id: usize,                   // ID of the first atom
    atom2_id: usize,                   // ID of the central atom
    atom3_id: usize,                   // ID of the third atom
    equilibrium_angle: f64,            // Equilibrium bond angle (θ0) in degrees
    force_constant: f64,               // Force constant (k_θ) in kcal/mol·rad²
}

#[derive(Debug, Clone)]
struct Torsion {
    atom1_id: usize,                   // ID of the first atom
    atom2_id: usize,                   // ID of the second atom (bonded to atom1)
    atom3_id: usize,                   // ID of the third atom (bonded to atom2)
    atom4_id: usize,                   // ID of the fourth atom (bonded to atom3)
    barrier_height: f64,               // Barrier height (V_n) in kcal/mol
    periodicity: u32,                  // Periodicity (n)
    phase: f64,                        // Phase angle (γ) in degrees
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

