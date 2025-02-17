/*
 * Digital Nets Conformational Sampling (DNCS)
 * Copyright [2024] [Abraham Rebairo J., Satheeshkumar S, Sam Paul D., Stephen A. and FiveAngstromLLP]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#![allow(dead_code)]
use clap::ValueEnum;
use serde::Deserialize;
use std::fmt::Debug;
use std::str::FromStr;
use std::sync::LazyLock;

pub static ALLAMINOMOLS_1: LazyLock<Polymer<Atom>> = LazyLock::new(|| {
    serde_json::from_str(include_str!("../library/jsons/ALLAMINOMOLS_1.json")).unwrap()
});
pub static ALLAMINOMOLS_2: LazyLock<Polymer<Atom>> = LazyLock::new(|| {
    serde_json::from_str(include_str!("../library/jsons/ALLAMINOMOLS_2.json")).unwrap()
});
pub static ALLCONN: LazyLock<Polymer<Bond>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../library/jsons/ALLCONN.json")).unwrap());
pub static DIHEDS: LazyLock<Polymer<Diheds>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../library/jsons/DIHEDS.json")).unwrap());
pub static ENERGYPARAM: LazyLock<EnergyParam> = LazyLock::new(|| {
    serde_json::from_str(include_str!("../library/jsons/ENERGYPARAM.json")).unwrap()
});
pub static VAR: LazyLock<Polymer<Diheds>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../library/jsons/VAR.json")).unwrap());

#[derive(Debug, Clone, ValueEnum)]
pub enum FF {
    Amber03,
    Amber10,
    Amber96,
    Amber99SB,
    AmberFB15,
}

impl FromStr for FF {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim_end_matches(".xml").to_lowercase().as_str() {
            "amber03" => Ok(FF::Amber03),
            "amber10" => Ok(FF::Amber10),
            "amber96" => Ok(FF::Amber96),
            "amber99-sb" => Ok(FF::Amber99SB),
            "amber-fb15" => Ok(FF::AmberFB15),
            _ => Err(format!("'{}' is not a valid forcefield", s)),
        }
    }
}

impl FF {
    pub fn init(&self) -> ForceField {
        match self {
            FF::Amber03 => {
                quick_xml::de::from_str(include_str!("../library/ForceFields/amber03.xml")).unwrap()
            }
            FF::Amber10 => {
                quick_xml::de::from_str(include_str!("../library/ForceFields/amber10.xml")).unwrap()
            }
            FF::Amber96 => {
                quick_xml::de::from_str(include_str!("../library/ForceFields/amber96.xml")).unwrap()
            }
            FF::Amber99SB => {
                quick_xml::de::from_str(include_str!("../library/ForceFields/amber99sb.xml"))
                    .unwrap()
            }
            FF::AmberFB15 => {
                quick_xml::de::from_str(include_str!("../library/ForceFields/amberfb15.xml"))
                    .unwrap()
            }
        }
    }
}

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
    pub mass: f64,
    /// Atomic charge
    pub charge: f64,
    /// Atom type (for force field calculations)
    pub atomtype: Option<String>,
    /// Atom TypeId
    pub typeid: Option<String>,
    /// Lennard-Jones parameter sigma
    pub sigma: f64,
    /// Lennard-Jones parameter epsilon
    pub epsilon: f64,
    /// Atomic velocity (vx, vy, vz)
    pub velocity: [f64; 3],
    /// Force acting on the atom (fx, fy, fz)
    pub force: [f64; 3],
}

impl Atom {
    pub fn new(line: String) -> Self {
        Atom {
            record: line.get(0..5).unwrap_or("").trim().to_string(),
            serial: line.get(7..12).unwrap_or("").trim().parse().unwrap_or(0),
            name: line.get(12..17).unwrap_or("").trim().to_string(),
            residue: line.get(17..21).unwrap_or("").trim().to_string(),
            chain_id: line.get(21..23).unwrap_or("").trim().to_string(),
            sequence: line.get(23..27).unwrap_or("").trim().parse().unwrap_or(0),
            position: [
                line.get(31..39).unwrap_or("").trim().parse().unwrap_or(0.0),
                line.get(39..47).unwrap_or("").trim().parse().unwrap_or(0.0),
                line.get(47..54).unwrap_or("").trim().parse().unwrap_or(0.0),
            ],
            occupancy: line.get(55..60).unwrap_or("").trim().parse().unwrap_or(0.0),
            bfactor: line.get(60..66).unwrap_or("").trim().parse().unwrap_or(0.0),
            element: line.get(76..78).unwrap_or("").trim().to_string(),
            mass: line.get(78..87).unwrap_or("").trim().parse().unwrap_or(0.0),
            charge: line.get(87..97).unwrap_or("").trim().parse().unwrap_or(0.0),
            atomtype: Some(line.get(97..).unwrap_or("").trim().to_string()),
            typeid: None,
            epsilon: 0.0,
            sigma: 0.0,
            velocity: [0.0; 3],
            force: [0.0; 3],
        }
    }
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
pub struct Bond {
    pub a: String,
    pub b: String,
}

/// Diheds
#[derive(Deserialize, Clone)]
pub struct Diheds {
    pub d: String,
    pub a: String,
    pub b: String,
    pub c: String,
}

/// Energy
#[derive(Deserialize, Clone)]
pub struct Energy {
    pub atomtype: String,
    pub sigma: f64,
    pub epsilon: f64,
}

/// EnergyParam
#[derive(Deserialize, Clone)]
pub struct EnergyParam {
    pub seq: Vec<Energy>,
}

/// Monomer
#[derive(Deserialize, Clone)]
pub struct Monomer<T> {
    /// Three-letter code for the monomer type
    pub tcode: String,
    /// Single-letter code for the monomer type
    pub scode: String,
    /// Number of atoms in the monomer
    pub natom: usize,
    /// Vector of atoms that make up the monomer
    pub atoms: Vec<T>,
}

/// Polymer
#[derive(Deserialize)]
pub struct Polymer<T> {
    pub seq: Vec<Monomer<T>>,
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
                atom.sequence = resno + 1;
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
pub fn atoms_to_pdbstring(atoms: Vec<Atom>) -> String {
    let mut a: Vec<Atom> = atoms
        .iter()
        .filter(|atom| atom.residue != "OH")
        .cloned()
        .collect();
    if let Some(terminal) = atoms
        .iter()
        .find(|atom| atom.residue == "OH" && atom.name == "O1")
        .cloned()
    {
        let mut terminal = terminal;
        terminal.name = "OXT".to_string();
        terminal.residue = a.last().unwrap().residue.clone();
        a.push(terminal);
    }
    a.iter_mut().enumerate().for_each(|(i, atom)| {
        atom.serial = i;
    });
    a.remove(0);
    a.iter()
        .map(|atom| format!("{:?}", atom))
        .collect::<Vec<String>>()
        .join("\n")
}

/// PDB to Atom
pub fn pdb_to_atoms(pdb_string: &str) -> Vec<Atom> {
    std::fs::read_to_string(pdb_string)
        .expect("Failed to read PDB file")
        .lines()
        .filter(|line| line.starts_with("ATOM") || line.starts_with("HETATM"))
        .map(|line| Atom::new(line.to_string()))
        .map(|atom| {
            if atom.name == "OXT" {
                let mut newatom = atom.clone();
                newatom.name = "HO".to_string();
                if newatom.sequence % 2 == 1 {
                    newatom.position[0] += 0.250;
                    newatom.position[1] += 0.927;
                    newatom.position[2] += 0.013;
                } else {
                    newatom.position[0] += 0.935;
                    newatom.position[1] += -0.218;
                    newatom.position[2] += -0.003;
                }
                newatom
            } else {
                atom
            }
        })
        .map(|atom| {
            if atom.name.len() == 3
                && (atom.name.contains("HA")
                    || atom.name.contains("HB")
                    || atom.name.contains("HG")
                    || atom.name.contains("HD"))
            {
                let d = atom.name.chars().last().unwrap();
                if d.is_numeric() {
                    let mut newatom = atom;
                    let mdigit = 4 - d.to_digit(10).unwrap();
                    newatom.name = newatom.name.replace(d, &mdigit.to_string());
                    newatom
                } else {
                    atom
                }
            } else {
                atom
            }
        })
        .collect()
}

/// Atoms to Sequence
pub fn atoms_to_seq(atoms: Vec<Atom>) -> String {
    let mut seq = "".to_string();
    let mut n = 0;
    for atom in atoms.iter() {
        if atom.sequence != n {
            match atom.residue.as_str() {
                "ALA" => seq.push('A'),
                "ARG" => seq.push('R'),
                "ASN" => seq.push('N'),
                "ASP" => seq.push('D'),
                "CYS" => seq.push('C'),
                "GLU" => seq.push('E'),
                "GLN" => seq.push('Q'),
                "GLY" => seq.push('G'),
                "HIS" => seq.push('H'),
                "ILE" => seq.push('I'),
                "LEU" => seq.push('L'),
                "LYS" => seq.push('K'),
                "MET" => seq.push('M'),
                "PHE" => seq.push('F'),
                "PRO" => seq.push('P'),
                "SER" => seq.push('S'),
                "THR" => seq.push('T'),
                "TRP" => seq.push('W'),
                "TYR" => seq.push('Y'),
                "VAL" => seq.push('V'),
                "AIB" => seq.push('B'),
                _ => {}
            }
            n += 1;
        }
    }
    seq
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "PascalCase")]
pub struct ForceField {
    pub atom_types: AtomTypes,
    pub residues: Residues,
    pub harmonic_bond_force: HarmonicBondForce,
    pub harmonic_angle_force: HarmonicAngleForce,
    pub periodic_torsion_force: PeriodicTorsionForce,
    pub nonbonded_force: NonbondedForce,
}

#[derive(Debug, Deserialize, Clone)]
pub struct AtomTypes {
    #[serde(rename = "Type")]
    pub types: Vec<AtomType>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct AtomType {
    #[serde(rename = "@name")]
    pub name: usize,
    #[serde(rename = "@class")]
    pub class: String,
    #[serde(rename = "@element")]
    pub element: String,
    #[serde(rename = "@mass")]
    pub mass: f64,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "PascalCase")]
pub struct Residues {
    pub residue: Vec<Residue>,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "PascalCase")]
pub struct Residue {
    #[serde(rename = "@name")]
    pub name: String,
    pub atom: Vec<ResidueAtom>,
    pub bond: Option<Vec<ResidueBond>>,
    pub external_bond: Option<Vec<ExternalBond>>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ResidueAtom {
    #[serde(rename = "@name")]
    pub name: String,
    #[serde(rename = "@type")]
    pub atype: usize,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ResidueBond {
    #[serde(rename = "@from")]
    pub from: usize,
    #[serde(rename = "@to")]
    pub to: usize,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ExternalBond {
    #[serde(rename = "@from")]
    pub from: usize,
}

#[derive(Deserialize, Debug, Clone)]
pub struct HarmonicBondForce {
    #[serde(rename = "Bond")]
    pub bonds: Vec<HarmonicBond>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct HarmonicBond {
    #[serde(rename = "@class1")]
    pub class1: String,
    #[serde(rename = "@class2")]
    pub class2: String,
    #[serde(rename = "@length")]
    pub length: f64,
    #[serde(rename = "@k")]
    pub k: f64,
}

#[derive(Deserialize, Debug, Clone)]
pub struct HarmonicAngleForce {
    #[serde(rename = "Angle")]
    pub angles: Vec<HarmonicAngle>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct HarmonicAngle {
    #[serde(rename = "@class1")]
    pub class1: String,
    #[serde(rename = "@class2")]
    pub class2: String,
    #[serde(rename = "@class3")]
    pub class3: String,
    #[serde(rename = "@angle")]
    pub angle: f64,
    #[serde(rename = "@k")]
    pub k: f64,
}

#[derive(Deserialize, Debug, Clone)]
pub struct PeriodicTorsionForce {
    #[serde(rename = "Proper")]
    pub proper: Vec<TorsionEntry>,
    #[serde(rename = "Improper")]
    pub improper: Vec<TorsionEntry>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct TorsionEntry {
    #[serde(rename = "@class1")]
    pub class1: String,
    #[serde(rename = "@class2")]
    pub class2: String,
    #[serde(rename = "@class3")]
    pub class3: String,
    #[serde(rename = "@class4")]
    pub class4: String,
    #[serde(rename = "@periodicity1")]
    pub periodicity1: f64,
    #[serde(rename = "@phase1")]
    pub phase1: f64,
    #[serde(rename = "@k1")]
    pub k1: f64,
    #[serde(rename = "@periodicity2")]
    pub periodicity2: Option<f64>,
    #[serde(rename = "@phase2")]
    pub phase2: Option<f64>,
    #[serde(rename = "@k2")]
    pub k2: Option<f64>,
    #[serde(rename = "@periodicity3")]
    pub periodicity3: Option<f64>,
    #[serde(rename = "@phase3")]
    pub phase3: Option<f64>,
    #[serde(rename = "@k3")]
    pub k3: Option<f64>,
    #[serde(rename = "@periodicity4")]
    pub periodicity4: Option<f64>,
    #[serde(rename = "@phase4")]
    pub phase4: Option<f64>,
    #[serde(rename = "@k4")]
    pub k4: Option<f64>,
    #[serde(rename = "@periodicity5")]
    pub periodicity5: Option<f64>,
    #[serde(rename = "@phase5")]
    pub phase5: Option<f64>,
    #[serde(rename = "@k5")]
    pub k5: Option<f64>,
    #[serde(rename = "@periodicity6")]
    pub periodicity6: Option<f64>,
    #[serde(rename = "@phase6")]
    pub phase6: Option<f64>,
    #[serde(rename = "@k6")]
    pub k6: Option<f64>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct NonbondedForce {
    #[serde(rename = "@coulomb14scale")]
    pub coulomb14scale: f64,
    #[serde(rename = "@lj14scale")]
    pub lj14scale: f64,
    #[serde(rename = "Atom")]
    pub atoms: Vec<NonbondedAtom>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct NonbondedAtom {
    #[serde(rename = "@type")]
    pub atom_type: usize,
    #[serde(rename = "@charge")]
    pub charge: f64,
    #[serde(rename = "@sigma")]
    pub sigma: f64,
    #[serde(rename = "@epsilon")]
    pub epsilon: f64,
}
