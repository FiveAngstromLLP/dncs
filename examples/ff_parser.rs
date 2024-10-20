#![allow(dead_code)]
use serde::Deserialize;

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Open the XML file
    let file = std::fs::File::open("data/amber99sb.xml")?;
    let reader = std::io::BufReader::new(file);

    // Deserialize the XML into our ForceField structure
    let force_field: ForceField = quick_xml::de::from_reader(reader)?;

    println!("{:?}", force_field);

    Ok(())
}
