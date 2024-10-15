#![allow(dead_code, non_snake_case)]

use libdncs::forcefield::Amber;
use libdncs::sampling::{RotateAtDihedral, Sampler};
use libdncs::system::System;
use nalgebra::Vector3;
use pyo3::prelude::*;

#[pyfunction]
fn getPDB(seq: String, filename: String) {
    let polymer = System::new(&seq);
    polymer.to_pdb(&filename);
}

#[pyfunction]
fn pdb_to_angle(filename: String, include_side_chain: bool) -> String {
    let angle = RotateAtDihedral::from_pdb(&filename, include_side_chain);
    let angle_csv = angle
        .iter()
        .map(|&a| a.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    angle_csv
}

#[pyfunction]
fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let u1 = Vector3::new(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    let u2 = Vector3::new(c[0] - b[0], c[1] - b[1], c[2] - b[2]);
    let u3 = Vector3::new(d[0] - c[0], d[1] - c[1], d[2] - c[2]);
    let sinth = (u2.norm() * u1).dot(&u2.cross(&u3));
    let costh = u1.cross(&u2).dot(&u2.cross(&u3));
    if u1.cross(&u2).dot(&u2.cross(&u3).cross(&u2)) < 0.0 {
        -sinth.atan2(costh).to_degrees()
    } else {
        sinth.atan2(costh).to_degrees()
    }
}

#[pyclass]
pub struct Polymer {
    pub polymer: System,
}

#[pymethods]
impl Polymer {
    #[new]
    fn fromAminoSEQ(seq: String) -> PyResult<Polymer> {
        let mut system = System::new(&seq);
        system.init_parameters();
        Ok(Polymer { polymer: system })
    }

    fn getEnergy(&mut self) -> PyResult<f64> {
        let ff = Amber::new(self.polymer.clone());
        Ok(ff.energy())
    }

    fn toPDB(&self, filename: String) {
        self.polymer.to_pdb(&filename);
    }

    fn dihedral(&mut self, sidechain: bool, foldername: String) {
        self.polymer.get_dihedralatoms(sidechain);
        self.polymer.dihedral_log(&foldername)
    }
}

#[pyclass]
pub struct SobolSampler {
    pub sampler: Sampler,
}

#[pymethods]
impl SobolSampler {
    #[new]
    fn new(system: &Polymer, no_of_samples: usize, sidechain: bool) -> PyResult<SobolSampler> {
        let mut sample = Sampler::new(system.polymer.clone());
        sample.system.get_dihedralatoms(sidechain);
        sample.sample(no_of_samples);
        Ok(SobolSampler { sampler: sample })
    }

    fn conformational_sort(&mut self) {
        self.sampler.conformational_sort();
    }

    fn toPDB(&self, filename: String) {
        self.sampler.to_pdb(&filename);
    }

    fn toPDBFiles(&self, prefix: String) {
        self.sampler.to_pdbfiles(&prefix);
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn dncs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(getPDB, m)?)?;
    m.add_function(wrap_pyfunction!(pdb_to_angle, m)?)?;
    m.add_function(wrap_pyfunction!(dihedral_angle, m)?)?;
    m.add_class::<Polymer>()?;
    m.add_class::<SobolSampler>()?;
    Ok(())
}
