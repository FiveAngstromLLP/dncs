#![allow(dead_code)]
use serde::Deserialize;

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "PascalCase")]
pub struct GlycamForceField {
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
    pub name: String,
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
    #[serde(rename = "@override")]
    pub _override: usize,
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
    #[serde(rename = "@charge")]
    pub charge: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ResidueBond {
    #[serde(rename = "@atomName1")]
    pub from: String,
    #[serde(rename = "@atomName2")]
    pub to: String,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ExternalBond {
    #[serde(rename = "@atomName")]
    pub from: String,
}

#[derive(Deserialize, Debug, Clone)]
pub struct HarmonicBondForce {
    #[serde(rename = "Bond")]
    pub bonds: Vec<HarmonicBond>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct HarmonicBond {
    #[serde(rename = "@type1")]
    pub class1: String,
    #[serde(rename = "@type2")]
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
    #[serde(rename = "@type1")]
    pub class1: String,
    #[serde(rename = "@type2")]
    pub class2: String,
    #[serde(rename = "@type3")]
    pub class3: String,
    #[serde(rename = "@angle")]
    pub angle: f64,
    #[serde(rename = "@k")]
    pub k: f64,
}

#[derive(Deserialize, Debug, Clone)]
pub struct PeriodicTorsionForce {
    #[serde(rename = "@ordering")]
    pub ordering: String,
    #[serde(rename = "Proper")]
    pub proper: Vec<TorsionEntry>,
    #[serde(rename = "Improper")]
    pub improper: Vec<TorsionEntry>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct TorsionEntry {
    #[serde(rename = "@type1")]
    pub class1: String,
    #[serde(rename = "@type2")]
    pub class2: String,
    #[serde(rename = "@type3")]
    pub class3: String,
    #[serde(rename = "@type4")]
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
