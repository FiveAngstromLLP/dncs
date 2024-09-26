use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::io::Write;
use std::sync::LazyLock;

const ALLAMINOMOLS_1: LazyLock<Polymer<Atom>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/ALLAMINOMOLS_1.json")).unwrap());
const ALLAMINOMOLS_2: LazyLock<Polymer<Atom>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/ALLAMINOMOLS_2.json")).unwrap());

/// Atom
#[derive(Deserialize, Serialize, Clone, PartialEq)]
pub struct Atom {
    pub record: String,
    pub serial: usize,
    pub name: String,
    pub residue: String,
    pub sequence: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[allow(dead_code)]
impl Atom {
    /// Create Atom from pdb formated line
    pub(crate) fn new(line: &str) -> Self {
        Atom {
            record: line.get(0..5).unwrap_or("").trim().to_string(),
            serial: line.get(7..12).unwrap_or("").trim().parse().unwrap_or(0),
            name: line.get(12..17).unwrap_or("").trim().to_string(),
            residue: line.get(17..21).unwrap_or("").trim().to_string(),
            sequence: line.get(23..27).unwrap_or("").trim().parse().unwrap_or(0),
            x: line.get(31..39).unwrap_or("").trim().parse().unwrap_or(0.0),
            y: line.get(39..47).unwrap_or("").trim().parse().unwrap_or(0.0),
            z: line.get(47..54).unwrap_or("").trim().parse().unwrap_or(0.0),
        }
    }
}

impl Debug for Atom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:<6}{:>5} {:^4} {:^4} {:>4}    {:>8.3}{:>8.3}{:>8.3}",
            self.record,
            self.serial,
            match self.name.chars().next() {
                Some(c) if c.is_alphabetic() && self.name.len() == 3 => format!(" {}", self.name),
                Some(c) if c.is_numeric() && self.name.len() == 3 => format!("{} ", self.name),
                _ => self.name.to_string(),
            },
            self.residue,
            self.sequence,
            self.x,
            self.y,
            self.z,
        )
    }
}

/// Bond
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct Bond {
    pub(crate) a: String,
    pub(crate) b: String,
}

/// Type
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct Type {
    pub(crate) name: String,
    pub(crate) symb: String,
    pub(crate) atomtype: String,
    pub(crate) charge: f64,
    pub(crate) val1: f64,
    pub(crate) val2: u32,
}

/// Diheds
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct Diheds {
    pub(crate) d: String,
    pub(crate) a: String,
    pub(crate) b: String,
    pub(crate) c: String,
}

/// Energy
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct Energy {
    pub(crate) atomtype: String,
    pub(crate) sigma: f64,
    pub(crate) epsilon: f64,
}

/// EnergyParam
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct EnergyParam {
    pub(crate) seq: Vec<Energy>,
}

/// Monomer
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct Monomer<T> {
    pub(crate) tcode: String,
    pub(crate) scode: String,
    pub(crate) natom: usize,
    pub(crate) atoms: Vec<T>,
}

/// Polymer
#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct Polymer<T> {
    pub(crate) seq: Vec<Monomer<T>>,
}

/// Polygen
pub struct Polygen {
    pub atoms: Vec<Atom>, // atoms
}

impl Polygen {
    /// Create Polymer
    pub fn new(seq: &str) -> Self {
        let seq = seq.trim().to_uppercase();
        for s in seq.chars() {
            assert!(!" OJUXZ".contains(s), "Invalid Sequence");
        }
        let mut sequence = Vec::<Atom>::new();
        let mut dumm = (0.0, 0.0, 0.0);
        let mut atomno = 1;
        let mut resno = 1;

        let data: Vec<(Monomer<Atom>, Monomer<Atom>)> = ALLAMINOMOLS_1
            .seq
            .par_iter()
            .zip(ALLAMINOMOLS_2.seq.par_iter())
            .map(|(a, b)| (a.clone(), b.clone()))
            .collect();

        if let Some(m) = data.par_iter().find_any(|lib| lib.0.tcode == "NH") {
            if !seq.starts_with('P') {
                for mut atom in m.0.atoms.clone() {
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
                        atomno += 1
                    } else if atom.record == "DUMM" {
                        continue;
                    }
                }
            } else {
                for mut atom in m.1.atoms.clone() {
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
                        atomno += 1
                    } else if atom.record == "DUMM" {
                        continue;
                    }
                }
            }
        }

        for (r, s) in seq.char_indices() {
            if let Some(m) = data.iter().find(|lib| lib.0.scode == s.to_string()) {
                if r % 2 == 0 {
                    for mut atom in m.0.atoms.clone() {
                        if atom.record == "DUMM" {
                            dumm.0 += atom.x;
                            dumm.1 += atom.y;
                            dumm.2 += atom.z;
                            break;
                        } else {
                            atom.serial = atomno;
                            atom.sequence = r + 1;
                            atom.x += dumm.0;
                            atom.y += dumm.1;
                            atom.z += dumm.2;
                            sequence.push(atom);
                            atomno += 1;
                            resno = r;
                        }
                    }
                } else {
                    for mut atom in m.1.atoms.clone() {
                        if atom.record == "DUMM" {
                            dumm.0 += atom.x;
                            dumm.1 += atom.y;
                            dumm.2 += atom.z;
                            break;
                        } else {
                            atom.serial = atomno;
                            atom.sequence = r + 1;
                            atom.x += dumm.0;
                            atom.y += dumm.1;
                            atom.z += dumm.2;
                            sequence.push(atom);
                            atomno += 1;
                            resno = r;
                        }
                    }
                }
            }
        }

        if let Some(m) = data.iter().find(|lib| lib.0.tcode == "COH") {
            if seq.len() % 2 == 1 {
                for mut atom in m.0.atoms.clone() {
                    if atom.record == "DUMM" {
                        dumm.0 += atom.x;
                        dumm.1 += atom.y;
                        dumm.2 += atom.z;
                        break;
                    } else {
                        atom.serial = atomno;
                        atom.sequence = resno + 1;
                        atom.x += dumm.0;
                        atom.y += dumm.1;
                        atom.z += dumm.2;
                        sequence.push(atom);
                        atomno += 1;
                    }
                }
            } else {
                for mut atom in m.1.atoms.clone() {
                    if atom.record == "DUMM" {
                        dumm.0 += atom.x;
                        dumm.1 += atom.y;
                        dumm.2 += atom.z;
                        break;
                    } else {
                        atom.serial = atomno;
                        atom.sequence = resno + 1;
                        atom.x += dumm.0;
                        atom.y += dumm.1;
                        atom.z += dumm.2;
                        sequence.push(atom);
                        atomno += 1;
                    }
                }
            }
        }
        Self { atoms: sequence }
    }

    /// Write pdb file
    #[allow(dead_code)]
    pub fn to_pdb(&self, filename: &str) -> Result<(), std::io::Error> {
        let mut outfile = std::fs::File::create(filename)?;
        for atom in &self.atoms {
            writeln!(outfile, "{:?}", atom)?;
        }
        Ok(())
    }
}
