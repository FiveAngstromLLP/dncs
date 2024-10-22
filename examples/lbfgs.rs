extern crate libdncs;

use libdncs::forcefield::Amber;
use libdncs::parser::AMBER99SB;
use libdncs::sampling::{RotateAtDihedral, Sampler};
use libdncs::system::System;
use liblbfgs::lbfgs;

fn main() {
    let mut s = System::new("AAAA", (*AMBER99SB).clone());
    s.init_parameters();
    s.get_dihedral();
    println!("{:?}", s.dihedral);
    let ff = Amber::new(s.clone());
    println!("Energy: {}", ff.energy());
    let mut sample = Sampler::new(s);
    sample.sample(10);
    let mut m = Minimizer::new(sample);
    m.minimize()
}

pub struct Minimizer {
    pub sample: Sampler,
    pub minimized: Vec<System>,
}

impl Minimizer {
    /// New Minimizer
    pub fn new(sample: Sampler) -> Self {
        Self {
            sample,
            minimized: Vec::new(),
        }
    }

    /// Automatic Differentiation of Energy Function
    #[inline]
    fn diff<'a>(func: Box<dyn Fn(Vec<f64>) -> f64 + 'a>, angle: Vec<f64>) -> (f64, Vec<f64>) {
        let h = 1.0 / 32.0;
        let fx: f64 = (func)(angle.clone());
        println!("diff {}", fx);
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
        for (sys, angle) in s.iter().zip(a) {
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
                .with_max_iterations(10)
                .minimize(&mut theta, evaluate, |prgr| {
                    println!("iter: {:?}; Energy: {}", prgr.niter, prgr.fx);
                    println!("x: {:?}", prgr.x);
                    false
                });
            match prbs {
                Ok(p) => {
                    println!("Minimized theta = {:?}", theta);
                    // self.minimized.sample_for_angle(theta);
                    println!("Minimized Energy = {:?}", p.fx);
                }
                Err(e) => {
                    println!("{:?}", e);
                }
            }
        }
    }
}
