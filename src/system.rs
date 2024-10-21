#![allow(dead_code)]
use crate::parser::{self, Atom, ForceField, ALLCONN, DIHEDS, ENERGYPARAM, VAR};
use std::io::Write;

pub type Particles = Vec<Atom>;

/// Represents a system of particles
#[derive(Debug, Clone)]
pub struct System {
    pub seq: String,
    pub forcefield: ForceField,
    pub particles: Particles,
    pub dihedral: Vec<(Atom, Atom, Atom, Atom)>,
    pub nonbonded: Vec<Particles>,
    pub bonded1_4: Vec<Particles>,
    pub hydrogen: Vec<(Atom, Atom)>,
}

impl System {
    /// Creates a new System from a sequence string
    pub fn new(seq: &str, forcefield: ForceField) -> Self {
        let sequence = parser::generate(seq);
        let total = sequence.len();
        Self {
            seq: seq.to_string(),
            forcefield,
            particles: sequence,
            dihedral: Vec::new(),
            nonbonded: vec![Vec::new(); total],
            bonded1_4: vec![Vec::new(); total],
            hydrogen: Vec::new(),
        }
    }

    pub fn from_pdb(file: &str, forcefield: ForceField) -> Self {
        let atoms: Vec<Atom>;
        atoms = parser::pdb_to_atoms(file);

        let total = atoms.len();
        Self {
            seq: parser::atoms_to_seq(atoms.clone()),
            forcefield,
            particles: atoms,
            dihedral: Vec::new(),
            nonbonded: vec![Vec::new(); total],
            bonded1_4: vec![Vec::new(); total],
            hydrogen: Vec::new(),
        }
    }

    pub fn get_atomtype(&self, atom: &Atom) -> Option<parser::AtomType> {
        for residue in &self.forcefield.residues.residue {
            if residue.name == atom.residue {
                for atom_type in &residue.atom {
                    if atom_type.name == atom.name {
                        for forcefield_type in &self.forcefield.atom_types.types {
                            if forcefield_type.name == atom_type.atype {
                                return Some(forcefield_type.clone());
                            }
                        }
                    }
                }
            }
        }
        None
    }

    pub fn init_parameters(&mut self) {
        self.get_neighbours();
        self.get_energyparameters();
        self.get_hydrogen_bonded();
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
    pub fn get_dihedralatoms(&mut self, sidechain: bool) {
        let mut val: (Atom, Atom, Atom, Atom) = (
            Atom::new("".to_string()),
            Atom::new("".to_string()),
            Atom::new("".to_string()),
            Atom::new("".to_string()),
        );
        for (i, s) in self.seq.chars().enumerate() {
            if let Some(m) = DIHEDS.seq.iter().find(|f| f.scode == s.to_string()) {
                let n = if sidechain { m.natom } else { 2 };
                for d in m.atoms.iter().take(n) {
                    for atom in self.particles.iter() {
                        if atom.name == d.a && atom.sequence == i + 1 {
                            val.0 = atom.clone()
                        }
                        if atom.name == d.b && atom.sequence == i + 1 {
                            val.1 = atom.clone()
                        }
                        if atom.name == d.c && atom.sequence == i + 1 {
                            val.2 = atom.clone()
                        }
                        if atom.name == d.d && atom.sequence == i + 1 {
                            val.3 = atom.clone()
                        }
                        if atom.name == d.d && atom.name == "HO" {
                            val.3 = atom.clone()
                        }
                    }
                    self.dihedral.push(val.clone())
                }
            }
        }
    }

    pub fn get_dihedral_for_angle(&mut self, sidechain: bool) {
        let lastseq = self.particles.last().unwrap().sequence;
        let mut val: (Atom, Atom, Atom, Atom) = (
            Atom::new("".to_string()),
            Atom::new("".to_string()),
            Atom::new("".to_string()),
            Atom::new("".to_string()),
        );

        for atom in self.particles.iter() {
            if atom.name == "H" && atom.sequence == 1 {
                val.0 = atom.clone();
            } else if atom.name == "N" && atom.sequence == 1 {
                val.1 = atom.clone();
            } else if atom.name == "CA" && atom.sequence == 1 {
                val.2 = atom.clone();
            } else if atom.name == "C" && atom.sequence == 1 {
                val.3 = atom.clone();
            }
        }
        self.dihedral.push(val.clone());

        for (i, s) in self.seq.chars().enumerate() {
            if let Some(m) = VAR.seq.iter().find(|f| f.scode == s.to_string()) {
                let n = if sidechain { m.natom } else { 2 };
                for d in m.atoms.iter().take(n) {
                    let seqb = if d.a == "C" && d.b == "N" && d.c == "CA" && d.d == "C" {
                        i + 2
                    } else {
                        i + 1
                    };
                    let seqd = match (d.a.as_str(), d.b.as_str(), d.c.as_str(), d.d.as_str()) {
                        ("C", "N", "CA", "C") | ("N", "CA", "C", "N") => i + 2,
                        _ => i + 1,
                    };

                    let get_atom = |name: &str, seq: usize| {
                        self.particles
                            .iter()
                            .find(|i| i.name == name && i.sequence == seq)
                    };

                    let a = get_atom(&d.a, i + 1);
                    let b = get_atom(&d.b, seqb);
                    let c = get_atom(&d.c, seqb);
                    let d = get_atom(&d.d, seqd);

                    if let (Some(a), Some(b), Some(c), Some(d)) = (a, b, c, d) {
                        self.dihedral
                            .push((a.clone(), b.clone(), c.clone(), d.clone()));
                    }
                }
            }
        }

        for atom in self.particles.iter() {
            if atom.name == "N" && atom.sequence == lastseq {
                val.0 = atom.clone();
            } else if atom.name == "CA" && atom.sequence == lastseq {
                val.1 = atom.clone();
            } else if atom.name == "C" && atom.sequence == lastseq {
                val.2 = atom.clone();
            } else if atom.name == "O" && atom.sequence == lastseq {
                val.3 = atom.clone();
            }
        }
        self.dihedral.push(val);
    }

    /// Find hydrogen-bonded atom pairs
    fn get_hydrogen_bonded(&mut self) {
        for i in self.particles.iter().take(self.particles.len() - 1) {
            if i.name == "H" || i.name == "HN" {
                for j in self.particles.iter().skip(i.serial + 1) {
                    if i.sequence < i.sequence && j.name == "O" {
                        self.hydrogen.push((i.clone(), j.clone()))
                    }
                }
            } else if i.name == "O" {
                for j in self.particles.iter().skip(i.serial + 1) {
                    if let Some(k) = self.particles.iter().find(|a| a.serial == i.serial + 1) {
                        if k.sequence < j.sequence && j.name == "H" {
                            self.hydrogen.push((i.clone(), j.clone()))
                        }
                    }
                }
            }
        }
    }

    /// Get Neighbours
    pub fn get_neighbours(&mut self) {
        for (i, atom) in self.particles.iter().enumerate() {
            let mut neighbor = Neighbor::new(self.particles.clone(), atom.clone());
            neighbor.get_neighbours();
            self.nonbonded[i] = neighbor.nonbonded;
            self.bonded1_4[i] = neighbor.bonded1_4;
        }
    }

    /// Writes the dihedral atoms
    pub fn dihedral_log(&self, foldername: &str) {
        let filename = format!("{}/dihedral.log", foldername);
        let mut file = std::fs::File::create(filename).unwrap();
        for (a, b, c, d) in self.dihedral.iter() {
            let atoma = self
                .particles
                .iter()
                .find(|i| i.serial == a.serial)
                .unwrap();
            let atomb = self
                .particles
                .iter()
                .find(|i| i.serial == b.serial)
                .unwrap();
            let atomc = self
                .particles
                .iter()
                .find(|i| i.serial == c.serial)
                .unwrap();
            let atomd = self
                .particles
                .iter()
                .find(|i| i.serial == d.serial)
                .unwrap();
            writeln!(
                file,
                "{} {} {} | {} {} {} | {} {} {} | {} {} {}",
                atoma.name,
                atoma.residue,
                atoma.sequence,
                atomb.name,
                atomb.residue,
                atomb.sequence,
                atomc.name,
                atomc.residue,
                atomc.sequence,
                atomd.name,
                atomd.residue,
                atomd.sequence,
            )
            .unwrap()
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

struct Neighbor {
    atom: Atom,
    polymer: Particles,
    nonbonded: Particles,
    bonded1_4: Particles,
}

impl Neighbor {
    pub fn new(polymer: Particles, atom: Atom) -> Self {
        Self {
            atom,
            polymer,
            nonbonded: Vec::new(),
            bonded1_4: Vec::new(),
        }
    }

    /// Get Neighbours
    fn get_neighbours(&mut self) {
        self.get_nonbonded();
        self.get_bonded1_4();
    }

    /// First Bonded atoms
    fn first_bonded(&self) -> Vec<Atom> {
        let mut first = Vec::new();
        for conn in ALLCONN.seq.iter().filter(|f| f.tcode == self.atom.residue) {
            for bond in conn.atoms.iter().filter(|b| b.a == self.atom.name) {
                for jatom in self.polymer.iter().filter(|j| j.name == bond.b) {
                    let itemp = match (bond.a.as_str(), bond.b.as_str()) {
                        ("C", "N") => self.atom.sequence + 1,
                        _ => self.atom.sequence,
                    };
                    if jatom.sequence == itemp && !first.contains(jatom) {
                        first.push(jatom.clone())
                    }
                }
            }
        }
        first
    }

    /// Second Bonded atoms
    fn second_bonded(&self) -> (Vec<Atom>, Vec<Atom>) {
        let first = self.first_bonded();
        let mut second = Vec::new();
        for conn in ALLCONN.seq.iter().filter(|f| f.tcode == self.atom.residue) {
            for fatom in first.iter() {
                for bond in conn.atoms.iter().filter(|b| fatom.name == b.a) {
                    for jatom in self.polymer.iter().filter(|j| j.name == bond.b) {
                        let itemp = match (bond.a.as_str(), bond.b.as_str()) {
                            ("C", "N") => self.atom.sequence + 1,
                            ("N", _) => fatom.sequence,
                            _ => self.atom.sequence,
                        };
                        if jatom.sequence == itemp
                            && self.atom.serial != jatom.serial
                            && !second.contains(jatom)
                        {
                            second.push(jatom.clone())
                        }
                    }
                }
            }
        }
        (first, second)
    }

    /// Third Bonded atoms
    fn third_bonded(&self) -> (Vec<Atom>, Vec<Atom>, Vec<Atom>) {
        let (first, second) = self.second_bonded();
        let mut third = Vec::new();
        for satom in second.iter() {
            for conn in ALLCONN.seq.iter().filter(|f| f.tcode == satom.residue) {
                for bond in conn.atoms.iter().filter(|b| satom.name == b.a) {
                    for jatom in self.polymer.iter().filter(|j| j.name == bond.b) {
                        let itemp = match (bond.a.as_str(), bond.b.as_str()) {
                            ("C", "N") => self.atom.sequence + 1,
                            ("N", _) | ("CA", _) => satom.sequence,
                            _ => self.atom.sequence,
                        };
                        if jatom.sequence == itemp
                            && self.atom.serial != jatom.serial
                            && jatom.serial > self.atom.serial
                            && !third.contains(jatom)
                            && !first.contains(jatom)
                        {
                            third.push(jatom.clone());
                        }
                    }
                }
            }
        }
        (first, second, third)
    }

    /// Fourth Bonded atoms
    fn fourth_bonded(&self) -> Vec<Atom> {
        let (first, second, third) = self.third_bonded();
        self.polymer
            .iter()
            .skip(self.atom.serial)
            .filter(|a| !first.contains(a) && !second.contains(a) && !third.contains(a))
            .cloned()
            .collect()
    }

    fn formated_fourth(&self) -> Vec<Atom> {
        let fourth = self.fourth_bonded();
        let mut val = Vec::new();
        if let Some(f) = fourth.first() {
            val.push(f.clone());
        }
        for atom in fourth.windows(2) {
            if atom[0].serial + 1 != atom[1].serial {
                val.push(atom[0].clone());
                val.push(atom[1].clone());
            }
        }
        if let Some(l) = fourth.last() {
            val.push(l.clone());
        }
        val
    }

    /// Calculate 1-4_bonded terms for given atom
    fn get_bonded1_4(&mut self) {
        let (first, second) = self.second_bonded();
        let mut bonded = Vec::new();
        for satom in second {
            for conn in ALLCONN.seq.iter().filter(|f| f.tcode == satom.residue) {
                for bond in conn.atoms.iter().filter(|b| satom.name == b.a) {
                    let itemp = match (bond.a.as_str(), bond.b.as_str()) {
                        ("C", "N") => self.atom.sequence + 1,
                        ("N", _) | ("CA", _) => satom.sequence,
                        _ => self.atom.sequence,
                    };
                    for jatom in self.polymer.iter().filter(|a| a.name == bond.b.as_str()) {
                        match (
                            self.atom.residue.as_str(),
                            self.atom.name.as_str(),
                            jatom.name.as_str(),
                        ) {
                            ("TYR", "CB", "HD1")
                            | ("TYR", "CB", "HD2")
                            | ("TYR", "CB", "CE1")
                            | ("TYR", "CB", "CE2")
                            | ("TYR", "CG", "CZ")
                            | ("TYR", "CG", "HE1")
                            | ("TYR", "CG", "HE2")
                            | ("TYR", "CE1", "HE2")
                            | ("TYR", "CE1", "HH")
                            | ("TYR", "CE2", "HH")
                            | ("TYR", "CD1", _)
                            | ("TYR", "HD1", _)
                            | ("TYR", "CD2", _)
                            | ("TYR", "HD2", _)
                            | ("TYR", "HE1", _)
                            | ("TYR", "HE2", _)
                            | ("PHE", "CB", "HD1")
                            | ("PHE", "CB", "HD2")
                            | ("PHE", "CB", "CE1")
                            | ("PHE", "CB", "CE2")
                            | ("PHE", "CG", "CZ")
                            | ("PHE", "CG", "HE1")
                            | ("PHE", "CG", "HE2")
                            | ("PHE", "CD1", _)
                            | ("PHE", "HD1", _)
                            | ("PHE", "CD2", _)
                            | ("PHE", "HD2", _)
                            | ("PHE", "CE1", _)
                            | ("PHE", "CE2", _)
                            | ("PHE", "HE1", _)
                            | ("PHE", "HE2", _)
                            | ("PHE", "CZ", _)
                            | ("PHE", "HZ", _)
                            | ("TRP", "CB", "HD1")
                            | ("TRP", "CB", "NE1")
                            | ("TRP", "CB", "CE2")
                            | ("TRP", "CB", "CE3")
                            | ("TRP", "CG", "CE2")
                            | ("TRP", "CG", "HE1")
                            | ("TRP", "CG", "CZ2")
                            | ("TRP", "CG", "HE3")
                            | ("TRP", "CG", "NE1")
                            | ("TRP", "CG", "CZ3")
                            | ("TRP", "CG", "HZ3")
                            | ("TRP", "CD1", _)
                            | ("TRP", "HD1", _)
                            | ("TRP", "CD2", _)
                            | ("TRP", "NE1", _)
                            | ("TRP", "HE1", _)
                            | ("TRP", "CE2", _)
                            | ("TRP", "CE3", _)
                            | ("TRP", "HE3", _)
                            | ("TRP", "CZ2", _)
                            | ("TRP", "HZ2", _)
                            | ("TRP", "CZ3", _)
                            | ("TRP", "HZ3", _)
                            | ("TRP", "CH2", _)
                            | ("TRP", "HH2", _)
                            | ("HIS", "CB", "CE1")
                            | ("HIS", "CB", "NE2")
                            | ("HIS", "CB", "HD2")
                            | ("HIS", "CG", "NE2")
                            | ("HIS", "CG", "HE1")
                            | ("HIS", "CG", "HE2")
                            | ("HIS", "CG", "CE1")
                            | ("HIS", "ND1", _)
                            | ("HIS", "CE1", _)
                            | ("HIS", "HE1", _)
                            | ("HIS", "NE2", _)
                            | ("HIS", "HE2", _)
                            | ("HIS", "CD2", _)
                            | ("HIS", "HD2", _)
                            | ("PRO", "C", _)
                            | ("PRO", "N", "1HB")
                            | ("PRO", "N", "2HB")
                            | ("PRO", "N", "CG")
                            | ("PRO", "N", "CB")
                            | ("PRO", "N", "1HG")
                            | ("PRO", "N", "2HG")
                            | ("PRO", "HA", "1HB")
                            | ("PRO", "HA", "2HB")
                            | ("PRO", "HA", "CG")
                            | ("PRO", "HA", "CB")
                            | ("PRO", "HA", "1HG")
                            | ("PRO", "HA", "2HG")
                            | ("PRO", "CB", "1HD")
                            | ("PRO", "CB", "2HD")
                            | ("PRO", "CB", "CD")
                            | ("PRO", "CA", _)
                            | ("PRO", "1HB", _)
                            | ("PRO", "2HB", _)
                            | ("PRO", "CG", _)
                            | ("PRO", "1HG", _)
                            | ("PRO", "2HG", _)
                            | ("PRO", "CD", _)
                            | ("PRO", "1HD", _)
                            | ("PRO", "2HD", _) => {
                                if self.atom.sequence == jatom.sequence {
                                    continue;
                                }
                            }
                            (_, "CA", "CA") | (_, "CA", "H") | (_, "O", "CA") | (_, "O", "H") => {
                                if self.atom.sequence + 1 == jatom.sequence {
                                    continue;
                                }
                            }
                            _ => {
                                if jatom.sequence == itemp
                                    && self.atom.serial != jatom.serial
                                    && jatom.serial > self.atom.serial
                                    && !bonded.contains(jatom)
                                    && !first.contains(jatom)
                                {
                                    bonded.push(jatom.clone());
                                }
                            }
                        }
                    }
                }
            }
        }
        self.bonded1_4 = bonded;
    }

    /// Calculate non-bonded terms for given atom
    fn get_nonbonded(&mut self) {
        let mut nonbonded = Vec::new();
        let mut skip_update = false;
        'fourth: for jatom in self.formated_fourth() {
            if self.atom.sequence == jatom.sequence {
                match (
                    self.atom.residue.as_str(),
                    self.atom.name.as_str(),
                    jatom.name.as_str(),
                ) {
                    ("TYR", "CB", "HE1")
                    | ("TYR", "CB", "HE2")
                    | ("TYR", "CB", "HH")
                    | ("TYR", "CG", "OH")
                    | ("TYR", "CG", "HH")
                    | ("PHE", "CB", "HE1")
                    | ("PHE", "CB", "HE2")
                    | ("PHE", "CB", "HZ")
                    | ("PHE", "CG", "HZ")
                    | ("TRP", "CB", "HE1")
                    | ("TRP", "CB", "HE3")
                    | ("TRP", "CB", "HH2")
                    | ("TRP", "CG", "HZ2")
                    | ("TRP", "CG", "HZ3")
                    | ("TRP", "CG", "HH2")
                    | ("HIS", "CB", "HE1")
                    | ("HIS", "CB", "HE2")
                    | ("PRO", "N", "1HG")
                    | ("PRO", "N", "2HG")
                    | ("PRO", "CA", "1HG")
                    | ("PRO", "CA", "2HG")
                    | ("PRO", "CB", "1HG")
                    | ("PRO", "CB", "2HG")
                    | ("PRO", "1HB", "2HG")
                    | ("PRO", "1HB", "1HG")
                    | ("PRO", "1HB", "2HD")
                    | ("PRO", "1HB", "1HD")
                    | ("PRO", "2HB", "2HG")
                    | ("PRO", "2HB", "1HG")
                    | ("PRO", "2HB", "2HD")
                    | ("PRO", "2HB", "1HD")
                    | ("PRO", "CG", "2HG")
                    | ("PRO", "CG", "1HG")
                    | ("PRO", "1HG", "2HG") => continue 'fourth,
                    ("TYR", "CD1", _)
                    | ("TYR", "HD1", _)
                    | ("TYR", "CD2", _)
                    | ("TYR", "HD2", _)
                    | ("TYR", "CE1", _)
                    | ("TYR", "HE1", _)
                    | ("TYR", "CE2", _)
                    | ("TYR", "HE2", _)
                    | ("TYR", "CZ", _)
                    | ("TYR", "OH", _)
                    | ("TYR", "HH", _)
                    | ("PHE", "CD1", _)
                    | ("PHE", "HD1", _)
                    | ("PHE", "CD2", _)
                    | ("PHE", "HD2", _)
                    | ("PHE", "CE1", _)
                    | ("PHE", "HE1", _)
                    | ("PHE", "CE2", _)
                    | ("PHE", "HE2", _)
                    | ("PHE", "CZ", _)
                    | ("PHE", "HZ", _)
                    | ("TRP", "CD1", _)
                    | ("TRP", "HD1", _)
                    | ("TRP", "CD2", _)
                    | ("TRP", "NE1", _)
                    | ("TRP", "HE1", _)
                    | ("TRP", "CE2", _)
                    | ("TRP", "CE3", _)
                    | ("TRP", "HE3", _)
                    | ("TRP", "CZ2", _)
                    | ("TRP", "HZ2", _)
                    | ("TRP", "CZ3", _)
                    | ("TRP", "HZ3", _)
                    | ("TRP", "CH2", _)
                    | ("TRP", "HH2", _)
                    | ("HIS", "HD2", _)
                    | ("HIS", "ND1", _)
                    | ("HIS", "CD2", _)
                    | ("HIS", "CE1", _)
                    | ("HIS", "HE1", _)
                    | ("HIS", "NE2", _)
                    | ("HIS", "HE2", _) => {
                        for l in self
                            .polymer
                            .iter()
                            .skip(self.atom.serial + 1)
                            .filter(|f| f.name == "C")
                        {
                            nonbonded.push(l.clone());
                            if let Some(last) = self.polymer.last() {
                                nonbonded.push(last.clone());
                            }
                            skip_update = true;
                        }
                    }
                    _ => nonbonded.push(jatom.clone()),
                }
            } else if self.atom.sequence + 1 == jatom.sequence && jatom.residue == "PRO" {
                match self.atom.name.as_str() {
                    "C" | "O" => {
                        for l in self
                            .polymer
                            .iter()
                            .skip(self.atom.serial + 1)
                            .filter(|f| f.name == "C")
                        {
                            nonbonded.push(l.clone());
                            if let Some(last) = self.polymer.last() {
                                nonbonded.push(last.clone());
                            }
                            skip_update = true;
                        }
                    }
                    _ => nonbonded.push(jatom.clone()),
                }
            } else {
                nonbonded.push(jatom.clone());
            }
        }
        if !skip_update {
            self.nonbonded = nonbonded
        }
    }
}
