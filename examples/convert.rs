use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, io::Write};

fn main() {
    convert()
}

fn convert() {
    // let data1: Polymer<Atom> = Polymer::from_lib("examples/files/ALLAMINOMOLS_1.lib");
    // data1.to_lib("examples/files/ALLAMINOMOLS_1.lib");
    // let data2: Polymer<Atom> = Polymer::from_lib("examples/files/ALLAMINOMOLS_2.lib");
    // data2.to_lib("examples/files/ALLAMINOMOLS_2.lib");
    // let data3: Polymer<Bond> = Polymer::from_lib("examples/files/ALLCONN.lib");
    // data3.to_lib("examples/files/ALLCONN.lib");
    // let data4: Polymer<Diheds> = Polymer::from_lib("examples/files/DIHEDS.lib");
    // data4.to_lib("examples/files/DIHEDS.lib");
    // let data5 = EnergyParam::from_lib("examples/files/ENERGYPARAM.lib");
    // data5.to_lib("examples/files/ENERGYPARAM.lib");
    // let data6: Polymer<Diheds> = Polymer::from_lib("examples/files/VAR.lib");
    // data6.to_lib("examples/files/VAR.lib");
    let data1: Polymer<Atom> = Polymer::from_lib("examples/files/ALLAMINOMOLS_1.lib");
    data1.to_json("data/ALLAMINOMOLS_1.json");
    let data2: Polymer<Atom> = Polymer::from_lib("examples/files/ALLAMINOMOLS_2.lib");
    data2.to_json("data/ALLAMINOMOLS_2.json");
    let data3: Polymer<Bond> = Polymer::from_lib("examples/files/ALLCONN.lib");
    data3.to_json("data/ALLCONN.json");
    let data4: Polymer<Diheds> = Polymer::from_lib("examples/files/DIHEDS.lib");
    data4.to_json("data/DIHEDS.json");
    let data5 = EnergyParam::from_lib("examples/files/ENERGYPARAM.lib");
    data5.to_json("data/ENERGYPARAM.json");
    let data6: Polymer<Diheds> = Polymer::from_lib("examples/files/VAR.lib");
    data6.to_json("data/VAR.json");
}

#[allow(dead_code)]
trait Convert {
    fn new(line: &str) -> Self;
    fn to_string(&self) -> String;
}

#[derive(Deserialize, Serialize, Clone)]
struct Atom {
    record: String,
    serial: usize,
    name: String,
    residue: String,
    chain_id: String,
    sequence: usize,
    position: [f64; 3],
    occupancy: f64,
    bfactor: f64,
    element: String,
    mass: f64,
    charge: f64,
    atomtype: String,
    sigma: f64,
    epsilon: f64,
    velocity: [f64; 3],
    force: [f64; 3],
}

impl Convert for Atom {
    fn new(line: &str) -> Self {
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
    fn to_string(&self) -> String {
        let name = convention(self.name.to_string());
        format!(
            "{:<6}{:>5} {:^4} {:^4} {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            self.record,
            self.serial,
            match self.name.len() {
                3 => format!(" {}", name),
                _ => name,
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

fn convention(name: String) -> String {
    if name.len() == 3
        && (name.contains("HA")
            || name.contains("HB")
            || name.contains("HG")
            || name.contains("HD"))
    {
        print!("{} ", name);
        let d = name.chars().last().unwrap();
        if d.is_numeric() {
            let mdigit = 4 - d.to_digit(10).unwrap();
            println!("{};", name.replace(d, &mdigit.to_string()));
            return name.replace(d, &mdigit.to_string());
        } else {
            return name;
        }
    } else {
        return name;
    }
}

impl Debug for Atom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = convention(self.name.to_string());
        write!(
            f,
            "{:<6}{:>5} {:^4} {:^4} {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            self.record,
            self.serial,
            match self.name.len() {
                3 => format!(" {}", name),
                _ => name,
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

#[derive(Deserialize, Serialize, Clone)]
struct Bond {
    a: String,
    b: String,
}

impl Convert for Bond {
    fn new(line: &str) -> Self {
        Bond {
            a: line.get(8..13).unwrap_or("").trim().to_string(),
            b: line.get(13..).unwrap_or("").trim().to_string(),
        }
    }
    fn to_string(&self) -> String {
        format!(" bond   {:^4} {:^4}", self.a, self.b)
    }
}

impl Debug for Bond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            " bond   {:^4} {:^4}",
            self.a.to_string(),
            self.b.to_string()
        )
    }
}

#[derive(Deserialize, Serialize, Clone)]
struct Energy {
    atomtype: String,
    sigma: f64,
    epsilon: f64,
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

impl Debug for Energy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "  {:<9}{:<9.3}{:<9.3}",
            self.atomtype, self.sigma, self.epsilon
        )
    }
}

#[derive(Deserialize, Serialize, Clone)]
struct Diheds {
    a: String,
    b: String,
    c: String,
    d: String,
}

impl Convert for Diheds {
    fn new(line: &str) -> Self {
        let var: Vec<String> = line
            .split(" ")
            .filter(|s| !s.is_empty())
            .map(String::from)
            .collect();
        Diheds {
            a: var[0].to_string(),
            b: var[1].to_string(),
            c: var[2].to_string(),
            d: var[3].to_string(),
        }
    }

    fn to_string(&self) -> String {
        let formt = |s: &str| match s.chars().next() {
            Some(c) if c.is_alphabetic() && s.len() == 3 => format!("{}", s.to_string()),
            Some(c) if c.is_numeric() && s.len() == 3 => format!("{}", s.to_string()),
            _ => s.to_string(),
        };
        format!(
            "        {:^6}{:^6}{:^6}{:^6}",
            formt(&self.a),
            formt(&self.b),
            formt(&self.c),
            formt(&self.d)
        )
    }
}

impl Debug for Diheds {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "        {:^6}{:^6}{:^6}{:^6}",
            self.a.to_string(),
            self.b.to_string(),
            self.c.to_string(),
            self.d.to_string()
        )
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct Monomer<T> {
    tcode: String,
    scode: String,
    natom: usize,
    atoms: Vec<T>,
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
struct EnergyParam {
    seq: Vec<Energy>,
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
