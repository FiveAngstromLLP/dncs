#![allow(dead_code, unused_imports)]
use crate::forcefield::Amber;
use crate::sampling::{RotateAtDihedral, Sampler};
use crate::system::System;
use liblbfgs::lbfgs;
use std::io::Write;

pub struct Minimizer {
    pub sample: Sampler,
    pub minimized: Vec<System>,
    pub energy: Vec<f64>,
    pub angles: Vec<Vec<f64>>,
}

impl Minimizer {
    /// New Minimizer
    pub fn new(sample: Sampler) -> Self {
        let len = sample.angles.len();
        Self {
            sample: sample.clone(),
            minimized: vec![sample.system.clone(); len],
            energy: vec![0.0; len],
            angles: vec![vec![0.0]; len],
        }
    }

    /// Automatic Differentiation of Energy Function
    #[inline]
    fn diff<'a>(func: Box<dyn Fn(Vec<f64>) -> f64 + 'a>, angle: Vec<f64>) -> (f64, Vec<f64>) {
        let h = 1.0 / 32.0;
        let fx: f64 = (func)(angle.clone());
        let dfx = angle
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let mut dx = angle.clone();
                dx[i] += h;
                ((func)(dx) - fx) / h
            })
            .collect();
        (fx, dfx)
    }

    /// Energy minimizer
    pub fn minimize(&mut self) {
        let s = self.sample.sample.clone();
        let a = self.sample.angles.clone();
        for (i, (sys, angle)) in s.iter().zip(a).enumerate() {
            let evaluate = |x: &[f64], gx: &mut [f64]| {
                let val = Self::diff(
                    Box::new(|x: Vec<f64>| RotateAtDihedral::new(sys.clone()).rotated_energy(x)),
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
                    println!("Model {} Energy : {:?} KCal/Mol", i, p.fx);
                    let mut r = RotateAtDihedral::new(self.sample.system.clone());
                    r.rotate(theta.clone());
                    self.minimized[i] = r.system.clone();
                    self.angles[i] = theta;
                    self.energy[i] = p.fx;
                }
                Err(e) => {
                    println!("{:?}", e);
                }
            }
        }
    }

    pub fn write_sampled_angles(&self, filename: &str) {
        let mut file = std::fs::File::create(filename).unwrap();
        for (i, (angles, energy)) in self.angles.iter().zip(self.energy.iter()).enumerate() {
            let line = format!(
                "{}, {} KCal/mol, {}\n",
                i + 1,
                energy,
                angles
                    .iter()
                    .map(|&a| a.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            );
            file.write_all(line.as_bytes()).unwrap();
        }
    }

    pub fn to_pdb(&self, filename: &str) {
        let mut file = std::fs::File::create(filename).unwrap();
        let eng = Amber::new(self.sample.system.clone()).energy();
        let pdb = RotateAtDihedral::new(self.sample.system.clone()).to_pdbstring(0 + 1, eng);
        file.write_all(pdb.as_bytes()).unwrap();
        for (i, s) in self.minimized.iter().enumerate() {
            let pdb = RotateAtDihedral::new(s.clone()).to_pdbstring(i + 1, self.energy[i]);
            file.write_all(pdb.as_bytes()).unwrap();
        }
    }
}
