#![allow(dead_code, non_snake_case)]

use libdncs::forcefield::Amber;
use libdncs::parser::AMBER99SB;
use libdncs::sampling::{RotateAtDihedral, Sampler};
use libdncs::system::System;
use pyo3::prelude::*;

#[pyfunction]
fn getPDB(seq: String, filename: String) {
    let polymer = System::new(&seq, (*AMBER99SB).clone());
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

#[pyclass]
pub struct Polymer {
    pub polymer: System,
}

#[pymethods]
impl Polymer {
    #[new]
    fn fromAminoSEQ(seq: String) -> PyResult<Polymer> {
        let mut system = System::new(&seq, (*AMBER99SB).clone());
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

    fn write_angles(&self, filename: String) {
        self.sampler.write_sampled_angles(&filename)
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
    m.add_class::<Polymer>()?;
    m.add_class::<SobolSampler>()?;
    Ok(())
}
