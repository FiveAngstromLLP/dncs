use crate::parser::{self, DIHEDS, ENERGYPARAM};
use std::io::Write;

pub type Particles = Vec<parser::Atom>;

/// Represents a system of particles
#[derive(Debug, Clone)]
pub struct System {
    pub seq: String,
    pub particles: Particles,
    pub dihedral: Vec<(usize, usize, usize, usize)>,
    pub nonbonded: Particles,
    pub bonded1_4: Particles,
}

impl System {
    /// Creates a new System from a sequence string
    pub fn new(seq: &str) -> Self {
        Self {
            seq: seq.to_string(),
            particles: parser::generate(seq),
            dihedral: Vec::new(),
            nonbonded: Vec::new(),
            bonded1_4: Vec::new(),
        }
    }

    pub fn init_parameters(&mut self) {
        self.get_energyparameters();
        self.get_dihedralatoms();
    }

    /// Get energy parameters
    fn get_energyparameters(&mut self) {
        self.particles.iter_mut().for_each(|atom| {
            if atom.serial == 1 {
                atom.sigma = 1.0000_f64;
                atom.epsilon = 0.0200_f64;
            } else {
                let typ = atom.atomtype.to_string();
                if atom.sequence == 1 && atom.name == "N" {
                    if let Some(enp) = ENERGYPARAM.seq.iter().find(|f| f.atomtype == "NT") {
                        atom.sigma = enp.sigma;
                        atom.epsilon = enp.epsilon;
                    }
                } else {
                    if let Some(enp) = ENERGYPARAM.seq.iter().find(|f| f.atomtype == typ) {
                        atom.sigma = enp.sigma;
                        atom.epsilon = enp.epsilon;
                    }
                }
            }
        });
    }

    /// Get Dihedral Atoms
    fn get_dihedralatoms(&mut self) {
        for (i, s) in self.seq.chars().enumerate() {
            if let Some(m) = DIHEDS.seq.iter().find(|f| f.scode == s.to_string()) {
                for d in m.atoms.iter().take(m.natom) {
                    let mut val = (0, 0, 0, 0);
                    for atom in self.particles.iter() {
                        if atom.name == d.a && atom.sequence == i + 1 {
                            val.0 = atom.serial
                        }
                        if atom.name == d.b && atom.sequence == i + 1 {
                            val.1 = atom.serial
                        }
                        if atom.name == d.c && atom.sequence == i + 1 {
                            val.2 = atom.serial
                        }
                        if atom.name == d.d && atom.sequence == i + 1 {
                            val.3 = atom.serial
                        }
                        if atom.name == d.d && atom.name == "HO" {
                            val.3 = atom.serial
                        }
                    }
                    self.dihedral.push(val)
                }
            }
        }
    }

    /// Writes the system to a PDB file
    pub fn to_pdb(&self, filename: &str) {
        let val = parser::atoms_to_pdbstring(self.particles.clone());
        let mut file = std::fs::File::create(filename).expect("Failed to create file");
        file.write_all(val.as_bytes())
            .expect("Failed to write to file");
    }
}
