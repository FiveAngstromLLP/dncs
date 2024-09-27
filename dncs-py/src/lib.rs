#![allow(dead_code, non_snake_case)]

use libdncs::forcefield::Amber;
use libdncs::sampling::Sampler;
use libdncs::system::System;
use pyo3::prelude::*;

// #[pyclass]
// pub struct PyAtom {
//     /// Record type (e.g., "ATOM" or "HETATM")
//     pub record: String,
//     /// Atom serial number
//     pub serial: usize,
//     /// Atom name
//     pub name: String,
//     /// Residue name
//     pub residue: String,
//     /// Chain identifier
//     pub chain_id: String,
//     /// Residue sequence number
//     pub sequence: usize,
//     /// Atom position (x, y, z coordinates)
//     pub position: [f64; 3],
//     /// Occupancy
//     pub occupancy: f64,
//     /// Temperature factor (B-factor)
//     pub bfactor: f64,
//     /// Element symbol
//     pub element: String,
// }

// #[pyfunction]
// fn fromAminoSEQ(seq: String) -> Vec<PyAtom> {
//     let atoms = parser::generate(&seq);
//     let mut val: Vec<PyAtom> = Vec::new();
//     for iatom in atoms {
//         val.push(PyAtom {
//             record: iatom.record,
//             serial: iatom.serial,
//             name: iatom.name.clone(),
//             residue: iatom.residue.clone(),
//             chain_id: iatom.chain_id.clone(),
//             sequence: iatom.sequence,
//             position: iatom.position,
//             occupancy: iatom.occupancy,
//             bfactor: iatom.bfactor,
//             element: iatom.element.clone(),
//         });
//     }
//     val
// }

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
}

#[pyclass]
pub struct SobolSampler {
    pub sampler: Sampler,
}

#[pymethods]
impl SobolSampler {
    #[new]
    fn new(system: &Polymer, no_of_samples: usize) -> PyResult<SobolSampler> {
        let mut sample = Sampler::new(system.polymer.clone());
        sample.sample(no_of_samples);
        Ok(SobolSampler { sampler: sample })
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
