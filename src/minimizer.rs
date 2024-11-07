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
use crate::{Amber, RotateAtDihedral, Sampler, System};
use liblbfgs::lbfgs;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

pub struct Minimizer {
    pub sample: Sampler,
    pub minimized: Vec<Arc<System>>,
    pub energy: Vec<f64>,
    pub angles: Vec<Vec<f64>>,
}

impl Minimizer {
    /// New Minimizer
    pub fn new(sample: Sampler) -> Self {
        let len = sample.angles.len();
        Self {
            sample: sample.clone(),
            minimized: vec![Arc::clone(&sample.system); len],
            energy: vec![0.0; len],
            angles: vec![vec![0.0]; len],
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

    /// Energy minimizer
    pub fn minimize(&mut self) {
        let s = self.sample.sample.clone();
        let a = self.sample.angles.clone();

        for (i, (sys, angle)) in s.iter().zip(a).enumerate() {
            let sys = sys.clone();
            let angle = angle.clone();
            let system_clone = self.sample.system.clone();

            let sys = Arc::new(sys.clone());
            let evaluate = |x: &[f64], gx: &mut [f64]| {
                let sys = Arc::clone(&sys);
                let val = Self::diff(
                    Arc::new(move |x: Vec<f64>| {
                        RotateAtDihedral::new((*sys).clone()).rotated_energy(x)
                    }),
                    x.to_vec(),
                );
                gx.copy_from_slice(&val.1);
                Ok(val.0)
            };

            let mut theta: Vec<f64> = angle.clone();
            let prbs = lbfgs()
                .with_max_iterations(5)
                .minimize(&mut theta, evaluate, |_| false);

            match prbs {
                Ok(p) => {
                    println!("Minimizing Model {}", i + 1);
                    let mut r = RotateAtDihedral::new(system_clone);
                    r.rotate(theta.clone());
                    self.minimized[i] = r.system;
                    self.angles[i] = theta;
                    self.energy[i] = p.fx;
                }
                Err(e) => {
                    println!("{:?}", e);
                }
            }
        }
    }

    pub fn conformational_sort(&mut self) {
        // const KBT: f64 = 300.0 * 1.380649e-23 * 6.02214076e23 / 4184.0; // KCal/mol

        // Create indices and sort them based on energy values
        let mut indices: Vec<usize> = (0..self.energy.len()).collect();
        indices.sort_by(|&a, &b| {
            self.energy[a]
                .partial_cmp(&self.energy[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Create new sorted vectors using the indices
        let sorted_energy: Vec<f64> = indices.iter().map(|&i| self.energy[i]).collect();
        let sorted_angles: Vec<Vec<f64>> =
            indices.iter().map(|&i| self.angles[i].clone()).collect();
        let sorted_minimized: Vec<Arc<System>> = indices
            .iter()
            .map(|&i| Arc::clone(&self.minimized[i]))
            .collect();

        // Replace original vectors with sorted ones
        let minenergy = sorted_energy.iter().copied().fold(f64::INFINITY, f64::min);
        self.energy = sorted_energy.iter().map(|eng| eng - minenergy).collect();
        self.angles = sorted_angles;
        self.minimized = sorted_minimized;
    }

    pub fn write_angles(&self, filename: &str) -> std::io::Result<()> {
        let mut file = File::create(filename)?;
        for (i, (angles, energy)) in self.angles.iter().zip(self.energy.iter()).enumerate() {
            let line = format!(
                "{}, {:<6.3}, {}\n",
                i + 1,
                energy,
                angles
                    .iter()
                    .map(|&a| a.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            );
            file.write_all(line.as_bytes())?;
        }
        Ok(())
    }

    pub fn to_pdb(&self, filename: &str) -> std::io::Result<()> {
        let mut file = File::create(filename)?;
        let eng = Amber::new(Arc::clone(&self.sample.system)).energy();
        let pdb = RotateAtDihedral::new(self.sample.system.clone()).to_pdbstring(1, eng);
        file.write_all(pdb.as_bytes())?;

        for (i, s) in self.minimized.iter().enumerate() {
            let pdb = RotateAtDihedral::new(Arc::clone(&s)).to_pdbstring(i + 1, self.energy[i]);
            file.write_all(pdb.as_bytes())?;
        }
        Ok(())
    }
}
