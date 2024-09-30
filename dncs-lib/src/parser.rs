#![allow(dead_code)]
use serde::Deserialize;
use std::fmt::Debug;
use std::sync::LazyLock;

pub(crate) static ALLAMINOMOLS_1: LazyLock<Polymer<Atom>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/ALLAMINOMOLS_1.json")).unwrap());
pub(crate) static ALLAMINOMOLS_2: LazyLock<Polymer<Atom>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/ALLAMINOMOLS_2.json")).unwrap());
pub(crate) static ALLCONN: LazyLock<Polymer<Bond>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/ALLCONN.json")).unwrap());
pub(crate) static DIHEDS: LazyLock<Polymer<Diheds>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/DIHEDS.json")).unwrap());
pub(crate) static ENERGYPARAM: LazyLock<EnergyParam> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/ENERGYPARAM.json")).unwrap());

#[derive(Deserialize, Clone, PartialEq)]
pub struct Atom {
    /// Record type (e.g., "ATOM" or "HETATM")
    pub record: String,
    /// Atom serial number
    pub serial: usize,
    /// Atom name
    pub name: String,
    /// Residue name
    pub residue: String,
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub sequence: usize,
    /// Atom position (x, y, z coordinates)
    pub position: [f64; 3],
    /// Occupancy
    pub occupancy: f64,
    /// Temperature factor (B-factor)
    pub bfactor: f64,
    /// Element symbol
    pub element: String,
    /// Atomic mass
    pub(crate) mass: f64,
    /// Atomic charge
    pub(crate) charge: f64,
    /// Atom type (for force field calculations)
    pub(crate) atomtype: String,
    /// Lennard-Jones parameter sigma
    pub(crate) sigma: f64,
    /// Lennard-Jones parameter epsilon
    pub(crate) epsilon: f64,
    /// Atomic velocity (vx, vy, vz)
    pub(crate) velocity: [f64; 3],
    /// Force acting on the atom (fx, fy, fz)
    pub(crate) force: [f64; 3],
}

impl Debug for Atom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:<6}{:>5} {:^4} {:^4} {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            self.record,
            self.serial,
            match self.name.len() {
                3 => format!(" {}", self.name),
                _ => self.name.to_string(),
            },
            self.residue,
            self.sequence,
            self.position[0],
            self.position[1],
            self.position[2],
            self.occupancy,
            self.bfactor,
            self.element
        )
    }
}

/// Bond
#[derive(Deserialize, Clone)]
pub(crate) struct Bond {
    pub(crate) a: String,
    pub(crate) b: String,
}

/// Diheds
#[derive(Deserialize, Clone)]
pub(crate) struct Diheds {
    pub(crate) d: String,
    pub(crate) a: String,
    pub(crate) b: String,
    pub(crate) c: String,
}

/// Energy
#[derive(Deserialize, Clone)]
pub(crate) struct Energy {
    pub(crate) atomtype: String,
    pub(crate) sigma: f64,
    pub(crate) epsilon: f64,
}

/// EnergyParam
#[derive(Deserialize, Clone)]
pub(crate) struct EnergyParam {
    pub(crate) seq: Vec<Energy>,
}

/// Monomer
#[derive(Deserialize, Clone)]
pub(crate) struct Monomer<T> {
    /// Three-letter code for the monomer type
    pub(crate) tcode: String,
    /// Single-letter code for the monomer type
    pub(crate) scode: String,
    /// Number of atoms in the monomer
    pub(crate) natom: usize,
    /// Vector of atoms that make up the monomer
    pub(crate) atoms: Vec<T>,
}

/// Polymer
#[derive(Deserialize)]
pub(crate) struct Polymer<T> {
    pub(crate) seq: Vec<Monomer<T>>,
}

/// Generates a polymer structure from a given amino acid sequence
pub fn generate(seq: &str) -> Vec<Atom> {
    let seq = seq.trim().to_uppercase();
    for s in seq.chars() {
        assert!(!" OJUXZ".contains(s), "Invalid Sequence");
    }

    let mut sequence = Vec::<Atom>::new();
    let mut dumm = [0.0; 3];
    let mut atomno = 1;
    let mut resno = 1;

    let data: Vec<(Monomer<Atom>, Monomer<Atom>)> = ALLAMINOMOLS_1
        .seq
        .iter()
        .zip(ALLAMINOMOLS_2.seq.iter())
        .map(|(a, b)| (a.clone(), b.clone()))
        .collect();

    // Handle N-terminal
    if let Some(m) = data.iter().find(|lib| lib.0.tcode == "NH") {
        let atoms = if !seq.starts_with('P') {
            &m.0.atoms
        } else {
            &m.1.atoms
        };
        for mut atom in atoms.clone() {
            if atom.record == "ATOM" {
                atom.serial = atomno;
                atom.sequence = 1;
                if let Some(r) = ALLAMINOMOLS_1
                    .seq
                    .iter()
                    .find(|f| f.scode.chars().next() == seq.chars().next())
                {
                    atom.residue = r.tcode.clone();
                }
                sequence.push(atom);
                atomno += 1;
            }
        }
    }

    // Process each amino acid in the sequence
    for (r, s) in seq.char_indices() {
        if let Some(m) = data.iter().find(|lib| lib.0.scode == s.to_string()) {
            let atoms = if r % 2 == 0 { &m.0.atoms } else { &m.1.atoms };
            for mut atom in atoms.clone() {
                if atom.record == "DUMM" {
                    dumm[0] += atom.position[0];
                    dumm[1] += atom.position[1];
                    dumm[2] += atom.position[2];
                    break;
                } else {
                    atom.serial = atomno;
                    atom.sequence = r + 1;
                    atom.position[0] += dumm[0];
                    atom.position[1] += dumm[1];
                    atom.position[2] += dumm[2];
                    sequence.push(atom);
                    atomno += 1;
                    resno = r;
                }
            }
        }
    }

    // Handle C-terminal
    if let Some(m) = data.iter().find(|lib| lib.0.tcode == "COH") {
        let atoms = if seq.len() % 2 == 1 {
            &m.0.atoms
        } else {
            &m.1.atoms
        };
        for mut atom in atoms.clone() {
            if atom.record == "DUMM" {
                dumm[0] += atom.position[0];
                dumm[1] += atom.position[1];
                dumm[2] += atom.position[2];
                break;
            } else {
                atom.serial = atomno;
                atom.sequence = resno + 2;
                atom.position[0] += dumm[0];
                atom.position[1] += dumm[1];
                atom.position[2] += dumm[2];
                sequence.push(atom);
                atomno += 1;
            }
        }
    }
    sequence
}

/// Write polymer structure to a PDB file
pub(crate) fn atoms_to_pdbstring(atoms: Vec<Atom>) -> String {
    let a: Vec<String> = atoms.iter().map(|atom| format!("{:?}", atom)).collect();
    // a.insert(a.len() - 2, "TER".to_string());
    a.join("\n")
}

// /// Atoms to Sequence
// pub(crate) fn atoms_to_seq(atoms: Vec<Atom>) -> String {
//     let mut seq = "".to_string();
//     let mut n = 0;
//     for atom in atoms.iter() {
//         if atom.sequence != n {
//             match atom.residue.as_str() {
//                 "ALA" => seq.push('A'),
//                 "ARG" => seq.push('R'),
//                 "ASN" => seq.push('N'),
//                 "ASP" => seq.push('D'),
//                 "CYS" => seq.push('C'),
//                 "GLU" => seq.push('E'),
//                 "GLN" => seq.push('Q'),
//                 "GLY" => seq.push('G'),
//                 "HIS" => seq.push('H'),
//                 "ILE" => seq.push('I'),
//                 "LEU" => seq.push('L'),
//                 "LYS" => seq.push('K'),
//                 "MET" => seq.push('M'),
//                 "PHE" => seq.push('F'),
//                 "PRO" => seq.push('P'),
//                 "SER" => seq.push('S'),
//                 "THR" => seq.push('T'),
//                 "TRP" => seq.push('W'),
//                 "TYR" => seq.push('Y'),
//                 "VAL" => seq.push('V'),
//                 "AIB" => seq.push('B'),
//                 _ => {}
//             }
//             n += 1;
//         }
//     }
//     seq
// }
