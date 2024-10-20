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
pub(crate) static VAR: LazyLock<Polymer<Diheds>> =
    LazyLock::new(|| serde_json::from_str(include_str!("../data/VAR.json")).unwrap());
pub(crate) static AMBER99SB: LazyLock<ForceField> =
    LazyLock::new(|| quick_xml::de::from_str(include_str!("../data/amber99sb.xml")).unwrap());

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
            atomtype: line.get(97..).unwrap_or("").trim().to_string(),
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
                // if newatom.sequence % 2 == 1 {
                //     newatom.position[0] += 0.250;
                //     newatom.position[1] += 0.927;
                //     newatom.position[2] += 0.013;
                // } else {
                //     newatom.position[0] += 0.935;
                //     newatom.position[1] += -0.218;
                //     newatom.position[2] += -0.003;
                // }
                return newatom;
            } else {
                return atom;
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
                    return newatom;
                } else {
                    return atom;
                }
            } else {
                return atom;
            }
        })
        .collect()
}

/// Atoms to Sequence
pub(crate) fn atoms_to_seq(atoms: Vec<Atom>) -> String {
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

#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase")]
pub(crate) struct ForceField {
    pub(crate) atom_types: AtomTypes,
    pub(crate) residues: Residues,
    pub(crate) harmonic_bond_force: HarmonicBondForce,
    pub(crate) harmonic_angle_force: HarmonicAngleForce,
    pub(crate) periodic_torsion_force: PeriodicTorsionForce,
    pub(crate) nonbonded_force: NonbondedForce,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase")]
pub(crate) struct Residues {
    pub(crate) residue: Vec<Residue>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase")]
pub(crate) struct Residue {
    #[serde(rename = "@name")]
    pub(crate) name: String,
    pub(crate) atom: Vec<ResidueAtom>,
    pub(crate) bond: Option<Vec<ResidueBond>>,
    pub(crate) external_bond: Option<Vec<ExternalBond>>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct ResidueAtom {
    #[serde(rename = "@name")]
    pub(crate) name: String,
    #[serde(rename = "@type")]
    pub(crate) atype: usize,
}

#[derive(Debug, Deserialize)]
pub(crate) struct ResidueBond {
    #[serde(rename = "@from")]
    pub(crate) from: usize,
    #[serde(rename = "@to")]
    pub(crate) to: usize,
}

#[derive(Debug, Deserialize)]
pub(crate) struct ExternalBond {
    #[serde(rename = "@from")]
    pub(crate) from: usize,
}

#[derive(Debug, Deserialize)]
pub(crate) struct AtomTypes {
    #[serde(rename = "Type")]
    pub(crate) types: Vec<AtomType>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct AtomType {
    #[serde(rename = "@name")]
    pub(crate) name: usize,
    #[serde(rename = "@class")]
    pub(crate) class: String,
    #[serde(rename = "@element")]
    pub(crate) element: String,
    #[serde(rename = "@mass")]
    pub(crate) mass: f64,
}

#[derive(Deserialize, Debug)]
pub(crate) struct HarmonicBondForce {
    #[serde(rename = "Bond")]
    pub(crate) bonds: Vec<HarmonicBond>,
}

#[derive(Deserialize, Debug)]
pub(crate) struct HarmonicBond {
    #[serde(rename = "@class1")]
    pub(crate) class1: String,
    #[serde(rename = "@class2")]
    pub(crate) class2: String,
    #[serde(rename = "@length")]
    pub(crate) length: f64,
    #[serde(rename = "@k")]
    pub(crate) k: f64,
}

#[derive(Deserialize, Debug)]
pub(crate) struct HarmonicAngleForce {
    #[serde(rename = "Angle")]
    pub(crate) angles: Vec<HarmonicAngle>,
}

#[derive(Deserialize, Debug)]
pub(crate) struct HarmonicAngle {
    #[serde(rename = "@class1")]
    pub(crate) class1: String,
    #[serde(rename = "@class2")]
    pub(crate) class2: String,
    #[serde(rename = "@class3")]
    pub(crate) class3: String,
    #[serde(rename = "@angle")]
    pub(crate) angle: f64,
    #[serde(rename = "@k")]
    pub(crate) k: f64,
}

#[derive(Deserialize, Debug)]
pub(crate) struct PeriodicTorsionForce {
    #[serde(rename = "Proper")]
    pub(crate) proper: Vec<ProperTorsion>,
    #[serde(rename = "Improper")]
    pub(crate) improper: Vec<ImproperTorsion>,
}

#[derive(Deserialize, Debug)]
pub(crate) struct ProperTorsion {
    #[serde(rename = "@class1")]
    pub(crate) class1: String,
    #[serde(rename = "@class2")]
    pub(crate) class2: String,
    #[serde(rename = "@class3")]
    pub(crate) class3: String,
    #[serde(rename = "@class4")]
    pub(crate) class4: String,
    #[serde(rename = "@periodicity1")]
    pub(crate) periodicity1: f64,
    #[serde(rename = "@phase1")]
    pub(crate) phase1: f64,
    #[serde(rename = "@k1")]
    pub(crate) k1: f64,
    #[serde(rename = "@periodicity2")]
    pub(crate) periodicity2: Option<f64>,
    #[serde(rename = "@phase2")]
    pub(crate) phase2: Option<f64>,
    #[serde(rename = "@k2")]
    pub(crate) k2: Option<f64>,
    #[serde(rename = "@periodicity3")]
    pub(crate) periodicity3: Option<f64>,
    #[serde(rename = "@phase3")]
    pub(crate) phase3: Option<f64>,
    #[serde(rename = "@k3")]
    pub(crate) k3: Option<f64>,
    #[serde(rename = "@periodicity4")]
    pub(crate) periodicity4: Option<f64>,
    #[serde(rename = "@phase4")]
    pub(crate) phase4: Option<f64>,
    #[serde(rename = "@k4")]
    pub(crate) k4: Option<f64>,
}

#[derive(Deserialize, Debug)]
pub(crate) struct ImproperTorsion {
    #[serde(rename = "@class1")]
    pub(crate) class1: String,
    #[serde(rename = "@class2")]
    pub(crate) class2: String,
    #[serde(rename = "@class3")]
    pub(crate) class3: String,
    #[serde(rename = "@class4")]
    pub(crate) class4: String,
    #[serde(rename = "@periodicity1")]
    pub(crate) periodicity1: f64,
    #[serde(rename = "@phase1")]
    pub(crate) phase1: f64,
    #[serde(rename = "@k1")]
    pub(crate) k1: f64,
}

#[derive(Deserialize, Debug)]
pub(crate) struct NonbondedForce {
    #[serde(rename = "@coulomb14scale")]
    pub(crate) coulomb14scale: f64,
    #[serde(rename = "@lj14scale")]
    pub(crate) lj14scale: f64,
    #[serde(rename = "Atom")]
    pub(crate) atoms: Vec<NonbondedAtom>,
}

#[derive(Deserialize, Debug)]
pub(crate) struct NonbondedAtom {
    #[serde(rename = "@type")]
    pub(crate) atom_type: usize,
    #[serde(rename = "@charge")]
    pub(crate) charge: f64,
    #[serde(rename = "@sigma")]
    pub(crate) sigma: f64,
    #[serde(rename = "@epsilon")]
    pub(crate) epsilon: f64,
}
