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

#![allow(dead_code, unused_imports)]
use crate::parser::{atoms_to_pdbstring, pdb_to_atoms, Atom, ForceField};
use crate::{forcefield, Amber, RotateAtDihedral, Sampler, System, FF};
use liblbfgs::lbfgs;
use std::fmt::format;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

pub struct Minimizer {
    pub samplefolder: String,
    pub folder: String,
    pub top_n_sample: usize,
    pub forcefield: ForceField,
    pub sample: Option<RotateAtDihedral>,
    pub energy: Vec<(usize, f64)>,
    pub angles: Vec<Vec<f64>>,
}

impl Minimizer {
    /// New Minimizer
    pub fn new(folder: String, forcefield: ForceField, top_n_sample: usize) -> Self {
        std::fs::create_dir_all(format!("{}/minimize", folder)).unwrap();
        Self {
            samplefolder: format!("{}/sample", folder),
            folder: format!("{}/minimize", folder),
            top_n_sample,
            forcefield,
            sample: None,
            energy: vec![(0, 0.0); top_n_sample],
            angles: vec![Vec::new(); top_n_sample],
        }
    }

    /// Automatic Differentiation of Energy Function
    #[inline]
    fn diff(
        func: Arc<dyn Fn(Vec<f64>) -> f64 + Send + Sync + 'static>,
        angle: Vec<f64>,
    ) -> (f64, Vec<f64>) {
        let h = 1.0 / 32.0;
        let fx: f64 = (func)(angle.clone());

        let mut dfx = Vec::with_capacity(angle.len());
        for (i, _) in angle.iter().enumerate() {
            let mut dx = angle.clone();
            dx[i] += h;
            let func = Arc::clone(&func);
            dfx.push(((func)(dx) - fx) / h);
        }

        (fx, dfx)
    }

    pub fn minimize_all(&mut self) {
        let files = std::fs::read_dir(&self.samplefolder).unwrap();
        let mut pdb_files: Vec<_> = files
            .filter_map(|entry| {
                let entry = entry.unwrap();
                let path = entry.path();
                if path.extension().map_or(false, |ext| ext == "pdb") {
                    let file_name = path.file_stem().unwrap().to_str().unwrap();
                    if file_name.starts_with("sample_") {
                        return Some(path);
                    }
                }
                None
            })
            .collect();

        pdb_files.sort_by_key(|path| {
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            let index_str = file_name.trim_start_matches("sample_");
            index_str.parse::<usize>().unwrap()
        });

        let top_pdb_files: Vec<String> = pdb_files
            .into_iter()
            .take(self.top_n_sample)
            .map(|path| path.to_string_lossy().into_owned())
            .collect();

        for (i, pdb_file) in top_pdb_files.iter().enumerate() {
            let mut sys = System::from_pdb(pdb_file, self.forcefield.clone());
            sys.init_parameters();
            self.minimize(Arc::new(sys), i);
        }

        // self.conformational_sort(300.0);
        // self.rename();
        self.write_angles();
    }

    /// Energy minimizer
    pub fn minimize(&mut self, system: Arc<System>, model: usize) {
        let system_clone = Arc::clone(&system);
        let evaluate = |x: &[f64], gx: &mut [f64]| {
            let sys = Arc::clone(&system_clone);
            let val = Self::diff(
                Arc::new(move |x: Vec<f64>| {
                    RotateAtDihedral::new(Arc::clone(&sys)).rotated_energy(x)
                }),
                x.to_vec(),
            );
            gx.copy_from_slice(&val.1);
            Ok(val.0)
        };

        let mut theta: Vec<f64> = RotateAtDihedral::new(Arc::clone(&system)).current_dihedral();
        let prbs = lbfgs()
            .with_max_iterations(5)
            .minimize(&mut theta, evaluate, |_| false);

        match prbs {
            Ok(p) => {
                println!("Minimizing Model {}", model + 1);
                let mut r = RotateAtDihedral::new(system_clone);
                r.rotate(theta.clone());
                let output_path = format!("{}/minimized_{:04}.pdb", self.folder, model);
                let mut output_file = File::create(output_path).unwrap();
                let pdb_string = atoms_to_pdbstring(r.rotated);
                output_file.write_all(pdb_string.as_bytes()).unwrap();
                self.angles[model] = theta;
                self.energy[model] = (model, p.fx);
            }
            Err(e) => {
                println!("{:?}", e);
            }
        }
    }

    pub fn conformational_sort(&mut self, temp: f64) {
        let _ = temp;
        // let kbt: f64 = temp * 1.380649e-23 * 6.02214076e23 / 4184.0; // KCal/mol
        // let weight: Vec<f64> = self.energy.iter().map(|e| (-e / kbt).exp()).collect();
        // let z: f64 = weight.iter().sum();
        // let normalized: Vec<f64> = weight.iter().map(|w| w / z).collect();

        // Create indices and sort them based on energy values

        self.energy.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let indices: Vec<usize> = self.energy.iter().map(|(i, _)| *i).collect();
        let angles: Vec<Vec<f64>> = (0..self.energy.len())
            .map(|i| self.angles[indices[i]].clone())
            .collect();
        self.angles = angles;
    }

    pub fn rename(&self) {
        // 1. Create a map of original indices to file paths:
        let file_map: std::collections::HashMap<usize, std::path::PathBuf> =
            std::fs::read_dir(&self.folder)
                .unwrap()
                .filter_map(|entry| {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    if path.extension().map_or(false, |ext| ext == "pdb") {
                        let file_name = path.file_stem().unwrap().to_str().unwrap();

                        if let Some(index_str) = file_name.split('_').last() {
                            if let Ok(index) = index_str.parse::<usize>() {
                                return Some((index, path));
                            }
                        }
                    }
                    None
                })
                .collect();

        // 2. Iterate over the sorted energy indices and rename files
        for (new_index, (original_index, _energy)) in self.energy.iter().enumerate() {
            if let Some(old_path) = file_map.get(original_index) {
                let new_path = format!("{}/minimze_{:04}.pdb", self.folder, new_index);
                if let Err(e) = std::fs::rename(old_path, &new_path) {
                    eprintln!("Failed to rename {:?} to {}: {}", old_path, new_path, e);
                }
            }
        }
    }

    pub fn write_angles(&self) {
        let mut file = std::fs::File::create(format!("{}/minimized.out", self.folder)).unwrap();
        for (i, (angles, energy)) in self.angles.iter().zip(self.energy.iter()).enumerate() {
            let line = format!(
                "{}, {:.3}, {}\n",
                i + 1,
                energy.1,
                angles
                    .iter()
                    .map(|&a| format!("{:.2}", a))
                    .collect::<Vec<String>>()
                    .join(", ")
            );
            file.write_all(line.as_bytes()).unwrap();
        }
    }
}
