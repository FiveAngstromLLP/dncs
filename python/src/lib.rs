/*
 * Digital Nets Conformational Sampling (DNCS)
 * Copyright [2024] [Abraham Rebairo J., Satheeshkumar S, Sam Paul D., Stephen A. and FiveAngstromLLP]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#![allow(dead_code, non_snake_case)]

use libdncs::*;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyfunction]
fn getPDB(seq: String, filename: String) {
    let polymer = System::new(&seq, FF::AMBERFB15.init());
    polymer.to_pdb(&filename);
}

#[pyfunction]
fn pdb_to_angle(filename: String) -> String {
    let angle = RotateAtDihedral::from_pdb(&filename);
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
    fn fromAminoSEQ(seq: String, forcefield: String) -> PyResult<Polymer> {
        // Validate forcefield
        let ff = match forcefield.as_str() {
            "amber03.xml" => FF::AMBER03,
            "amber10.xml" => FF::AMBER10,
            "amber96.xml" => FF::AMBER96,
            "amber99sb.xml" => FF::AMBER99SB,
            "amberfb15.xml" => FF::AMBERFB15,
            _ => {
                eprintln!(
                    "Unsupported forcefield: {}.
Must be one of the below:
  - amber03.xml
  - amber10.xml
  - amber96.xml
  - amber99sb.xml
  - amberfb15.xml
  ",
                    forcefield
                );
                std::process::exit(1);
            }
        };
        let mut system = System::new(&seq, ff.init());
        system.init_parameters();
        Ok(Polymer { polymer: system })
    }

    fn getEnergy(&mut self) -> PyResult<f64> {
        let ff = Amber::new(Arc::new(self.polymer.clone()));
        Ok(ff.energy())
    }

    fn toPDB(&self, filename: String) {
        self.polymer.to_pdb(&filename);
    }

    fn dihedral(&mut self, foldername: String) {
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
    fn new(
        system: &Polymer,
        no_of_samples: usize,
        grid: usize,
        temp: f64,
    ) -> PyResult<SobolSampler> {
        let mut sample = Sampler::new(Arc::new(system.polymer.clone()), grid);
        sample.sample(no_of_samples, temp);
        Ok(SobolSampler { sampler: sample })
    }

    fn write_angles(&self, filename: String) {
        self.sampler.write_angles(&filename)
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
