use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, io::Write};

fn main() {
    convert()
}

fn convert() {
    let data1: Polymer<Atom> = Polymer::from_lib("examples/files/ALLAMINOMOLS_1.lib");
    data1.to_json("data/ALLAMINOMOLS_1.json");
    let data2: Polymer<Atom> = Polymer::from_lib("examples/files/ALLAMINOMOLS_2.lib");
    data2.to_json("data/ALLAMINOMOLS_2.json");
    let data3: Polymer<Bond> = Polymer::from_lib("examples/files/ALLCONN.lib");
    data3.to_json("data/ALLCONN.json");
    let data4: Polymer<Type> = Polymer::from_lib("examples/files/ATOMTYPE.lib");
    data4.to_json("data/ATOMTYPE.json");
    let data5: Polymer<Diheds> = Polymer::from_lib("examples/files/DIHEDS.lib");
    data5.to_json("data/DIHEDS.json");
    let data6 = EnergyParam::from_lib("examples/files/ENERGYPARAM.lib");
    data6.to_json("data/ENERGYPARAM.json");
}

#[allow(dead_code)]
trait Convert {
    fn new(line: &str) -> Self;
    fn to_string(&self) -> String;
}

#[derive(Deserialize, Serialize, Clone)]
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

impl Convert for Atom {
    fn new(line: &str) -> Self {
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
    fn to_string(&self) -> String {
        format!(
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

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Bond {
    pub a: String,
    pub b: String,
}

impl Convert for Bond {
    fn new(line: &str) -> Self {
        Bond {
            a: line.get(8..12).unwrap_or("").trim().to_string(),
            b: line.get(12..).unwrap_or("").trim().to_string(),
        }
    }
    fn to_string(&self) -> String {
        format!(" bond   {:^4} {:^4}", self.a, self.b)
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Type {
    pub name: String,
    pub symb: String,
    pub atomtype: String,
    pub charge: f64,
    pub val1: f64,
    pub val2: u32,
}

impl Convert for Type {
    fn new(line: &str) -> Self {
        Type {
            name: line.get(0..5).unwrap_or("").trim().to_string(),
            symb: line.get(10..13).unwrap_or("").trim().to_string(),
            atomtype: line.get(13..17).unwrap_or("").trim().to_string(),
            charge: line.get(17..25).unwrap_or("").trim().parse().unwrap_or(0.0),
            val1: line.get(25..32).unwrap_or("").trim().parse().unwrap_or(0.0),
            val2: line.get(33..).unwrap_or("").trim().parse().unwrap_or(0),
        }
    }

    fn to_string(&self) -> String {
        format!(
            "{:^4}{:>8}  {:<3}{:>8.4}{:>8.3} {:>3}",
            match self.name.chars().next() {
                Some(c) if c.is_alphabetic() && self.name.len() == 3 => format!(" {}", self.name),
                Some(c) if c.is_numeric() && self.name.len() == 3 => format!("{} ", self.name),
                _ => self.name.to_string(),
            },
            self.symb,
            self.atomtype,
            self.charge,
            self.val1,
            self.val2
        )
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Energy {
    pub atomtype: String,
    pub sigma: f64,
    pub epsilon: f64,
}

impl Convert for Energy {
    fn new(line: &str) -> Self {
        let l: Vec<&str> = line.split_whitespace().collect();
        Energy {
            atomtype: l[0].to_string(),
            sigma: l[1].parse().unwrap_or(0.0),
            epsilon: l[2].parse().unwrap_or(0.0),
        }
    }

    fn to_string(&self) -> String {
        format!(
            "  {:<9}{:<9.3}{:<9.3}",
            self.atomtype, self.sigma, self.epsilon
        )
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Diheds {
    pub d: String,
    pub a: String,
    pub b: String,
    pub c: String,
}

impl Convert for Diheds {
    fn new(line: &str) -> Self {
        Diheds {
            a: line.get(9..13).unwrap_or("").trim().to_string(),
            b: line.get(13..16).unwrap_or("").trim().to_string(),
            c: line.get(16..20).unwrap_or("").trim().to_string(),
            d: line.get(20..).unwrap_or("").trim().to_string(),
        }
    }

    fn to_string(&self) -> String {
        let formt = |s: &str| match s.chars().next() {
            Some(c) if c.is_alphabetic() && s.len() == 3 => format!("{}", s.to_string()),
            Some(c) if c.is_numeric() && s.len() == 3 => format!("{}", s.to_string()),
            _ => s.to_string(),
        };
        format!(
            "        {:^4}{:^4}{:^4}{:^4}",
            formt(&self.a),
            formt(&self.b),
            formt(&self.c),
            formt(&self.d)
        )
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Var {
    pub a: String,
    pub b: String,
    pub c: String,
    pub d: String,
    pub t: String,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Convert for Var {
    fn new(line: &str) -> Self {
        Var {
            a: line.get(9..13).unwrap_or("").trim().to_string(),
            b: line.get(13..16).unwrap_or("").trim().to_string(),
            c: line.get(16..20).unwrap_or("").trim().to_string(),
            d: line.get(20..24).unwrap_or("").trim().to_string(),
            t: line.get(26..31).unwrap_or("").trim().to_string(),
            x: line.get(33..39).unwrap_or("").trim().parse().unwrap_or(0.0),
            y: line.get(39..46).unwrap_or("").trim().parse().unwrap_or(0.0),
            z: line.get(46..).unwrap_or("").trim().parse().unwrap_or(0.0),
        }
    }

    fn to_string(&self) -> String {
        let formt = |s: &str| match s.chars().next() {
            Some(c) if c.is_alphabetic() && s.len() == 3 => format!(" {}", s),
            Some(c) if c.is_numeric() && s.len() == 3 => format!("{} ", s),
            _ => s.to_string(),
        };
        format!(
            "        {:^4}{:^4}{:^4}{:^4}  {:<4}{:>7.2}{:>7.2}{:>7.2}",
            formt(&self.a),
            formt(&self.b),
            formt(&self.c),
            formt(&self.d),
            self.t,
            self.x,
            self.y,
            self.z
        )
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct SC {
    pub seq: String,
    pub serial: u32,
    pub val: Vec<f32>,
}

impl Convert for SC {
    fn new(line: &str) -> Self {
        SC {
            seq: line.get(0..4).unwrap_or("").trim().to_string(),
            serial: line.get(5..9).unwrap_or("").trim().parse().unwrap_or(0),
            val: line
                .split_whitespace()
                .skip(2)
                .map(|f| f.parse().unwrap_or(0.0))
                .collect(),
        }
    }
    fn to_string(&self) -> String {
        let mut val = "".to_string();
        for i in &self.val {
            val += format!("{:>8.1}", i).as_str()
        }
        format!("{:<5}{:>3}             {}", self.seq, self.serial, val)
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Monomer<T> {
    pub tcode: String,
    pub scode: String,
    pub natom: usize,
    pub atoms: Vec<T>,
}

impl<T> Convert for Monomer<T> {
    fn new(line: &str) -> Self {
        Monomer {
            tcode: line.get(0..4).unwrap_or("").trim().to_string(),
            scode: line.get(4..6).unwrap_or("").trim().to_string(),
            natom: line.get(6..).unwrap_or("").trim().parse().unwrap_or(0),
            atoms: Vec::new(),
        }
    }
    fn to_string(&self) -> String {
        "".to_string()
    }
}

#[allow(dead_code)]
trait FileHandling {
    fn from_lib(filepath: &str) -> Self;
    fn to_json(&self, filepath: &str);
    fn to_lib(&self, filepath: &str);
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct Polymer<T> {
    seq: Vec<Monomer<T>>,
}

impl<T: Convert + Serialize + Debug> FileHandling for Polymer<T> {
    fn from_lib(filepath: &str) -> Self {
        let data = std::fs::read_to_string(filepath).expect("File Not Found");
        let mut monomers = Vec::new();
        let lines: Vec<_> = data.lines().filter(|l| !l.is_empty()).collect();
        let ra = Regex::new(r"^[\sA-Z]+\s[A-Z\s]\s+\d+$").unwrap();
        let rb = Regex::new(r"^[\sA-Z]+\s[A-Z\s]\s+\d+\s+$").unwrap();
        let rc = Regex::new(r"^[\sA-Z]+\s[A-Z\s]\s+\d+\s+\d+$").unwrap();
        for line in lines {
            if ra.is_match(line) | rb.is_match(line) | rc.is_match(line) {
                monomers.push(Monomer::new(line))
            } else if let Some(last) = monomers.last_mut() {
                last.atoms.push(T::new(line));
            }
        }
        Polymer { seq: monomers }
    }

    fn to_json(&self, filepath: &str) {
        let data = serde_json::to_string_pretty(&self).unwrap();
        let mut file = std::fs::File::create(filepath).expect("File Not Created");
        file.write_all(data.as_bytes()).expect("Failed to write");
    }

    fn to_lib(&self, filepath: &str) {
        let mut outfile = std::fs::File::create(filepath).expect("File Not Created");
        for m in &self.seq {
            let l = format!("{:<4}{}{:>3}", m.tcode, m.scode, m.natom);
            writeln!(outfile, "{}", l).expect("Error Writing Monomer");
            for atom in &m.atoms {
                writeln!(outfile, "{:?}", atom).expect("Error Writing Atom")
            }
        }
    }
}

#[derive(Deserialize, Serialize, Clone)]
pub struct EnergyParam {
    pub seq: Vec<Energy>,
}

impl FileHandling for EnergyParam {
    fn from_lib(filepath: &str) -> Self {
        let data = std::fs::read_to_string(filepath).expect("File Not Found");
        let mut energyparam = Vec::new();
        let lines: Vec<_> = data.lines().filter(|l| !l.is_empty()).collect();
        for line in lines {
            energyparam.push(Energy::new(line))
        }
        EnergyParam { seq: energyparam }
    }

    fn to_json(&self, filepath: &str) {
        let data = serde_json::to_string_pretty(&self).unwrap();
        let mut file = std::fs::File::create(filepath).expect("File Not Created");
        file.write_all(data.as_bytes()).expect("Failed to write");
    }

    fn to_lib(&self, filepath: &str) {
        let mut outfile = std::fs::File::create(filepath).expect("File Not Created");
        for e in &self.seq {
            writeln!(outfile, "{}", e.to_string()).expect("Error Writing Atom")
        }
    }
}
