use crate::parser;
use std::io::Write;

pub type Particles = Vec<parser::Atom>;

/// Represents a system of particles
pub struct System {
    /// Vector of atoms in the system
    pub particles: Particles,
    pub samples: Vec<Particles>,
}

impl System {
    /// Creates a new System from a sequence string
    pub fn new(seq: &str) -> Self {
        Self {
            particles: parser::generate(seq),
            samples: Vec::new(),
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
