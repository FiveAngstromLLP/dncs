#![allow(dead_code, non_snake_case)]

use libdncs::forcefield::Amber;
use libdncs::sampling::Sampler;
use libdncs::system::System;
use pyo3::prelude::*;

#[pyfunction]
fn getPDB(seq: String, filename: String) {
    let polymer = System::new(&seq);
    polymer.to_pdb(&filename);
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

    fn dihedral(&self, foldername: String) {
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
    m.add_class::<Polymer>()?;
    m.add_class::<SobolSampler>()?;
    Ok(())
}
