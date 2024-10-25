#![allow(dead_code, unused_imports)]
use crate::{Amber, RotateAtDihedral, Sampler, System};
use liblbfgs::lbfgs;
use std::io::Write;
use std::sync::Arc;
use tokio::io::AsyncWriteExt;

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
    async fn diff(
        func: Arc<dyn Fn(Vec<f64>) -> f64 + Send + Sync + 'static>,
        angle: Vec<f64>,
    ) -> (f64, Vec<f64>) {
        let h = 1.0 / 32.0;
        let fx: f64 = (func)(angle.clone());

        let mut futures = Vec::new();
        for (i, _) in angle.iter().enumerate() {
            let mut dx = angle.clone();
            dx[i] += h;
            let func = Arc::clone(&func);
            futures.push(tokio::task::spawn_blocking(move || ((func)(dx) - fx) / h));
        }

        let mut dfx = Vec::with_capacity(angle.len());
        for future in futures {
            dfx.push(future.await.unwrap());
        }

        (fx, dfx)
    }

    /// Energy minimizer
    pub async fn minimize(&mut self) {
        let s = self.sample.sample.clone();
        let a = self.sample.angles.clone();

        let mut futures = Vec::new();
        for (i, (sys, angle)) in s.iter().zip(a).enumerate() {
            let sys = sys.clone();
            let angle = angle.clone();
            let system_clone = self.sample.system.clone();

            futures.push(tokio::task::spawn_blocking(move || {
                let sys = Arc::new(sys.clone());
                let evaluate = move |x: &[f64], gx: &mut [f64]| {
                    let sys = Arc::clone(&sys);
                    let val = futures::executor::block_on(Self::diff(
                        Arc::new(move |x: Vec<f64>| {
                            RotateAtDihedral::new((*sys).clone()).rotated_energy(x)
                        }),
                        x.to_vec(),
                    ));
                    gx.copy_from_slice(&val.1);
                    Ok(val.0)
                };

                let mut theta: Vec<f64> = angle.clone();
                let prbs = lbfgs()
                    .with_max_iterations(10)
                    .minimize(&mut theta, evaluate, |_| false);

                match prbs {
                    Ok(p) => {
                        println!("Model {} Energy : {:?} KCal/Mol", i + 1, p.fx);
                        let mut r = RotateAtDihedral::new(system_clone);
                        r.rotate(theta.clone());
                        Ok((r.system, theta, p.fx))
                    }
                    Err(e) => {
                        println!("{:?}", e);
                        Err(e)
                    }
                }
            }));
        }

        for (i, future) in futures.into_iter().enumerate() {
            if let Ok(Ok((system, theta, energy))) = future.await {
                self.minimized[i] = system;
                self.angles[i] = theta;
                self.energy[i] = energy;
            }
        }
    }

    pub async fn conformational_sort(&mut self) {
        const KBT: f64 = 300.0 * 1.380649e-23 * 6.02214076e23 / 4184.0; // KCal/mol
        let weight: Vec<f64> = self.energy.iter().map(|e| (-e / KBT).exp()).collect();
        let z: f64 = weight.iter().sum();
        let normalized: Vec<f64> = weight.iter().map(|w| w / z).collect();

        let energies = std::mem::take(&mut self.energy);
        let angles = std::mem::take(&mut self.angles);
        let minimized = std::mem::take(&mut self.minimized);

        let mut combined: Vec<(f64, Vec<f64>, System, f64)> = energies
            .into_iter()
            .zip(angles.into_iter())
            .zip(minimized.into_iter())
            .zip(normalized.into_iter())
            .map(|(((e, a), s), w)| (e, a, s, w))
            .collect();

        combined.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));

        self.energy = combined.iter().map(|(e, _, _, _)| *e).collect();
        self.angles = combined.iter().map(|(_, a, _, _)| a.clone()).collect();
        self.minimized = combined.into_iter().map(|(_, _, s, _)| s).collect();
    }

    pub async fn write_angles(&self, filename: &str) -> std::io::Result<()> {
        let mut file = tokio::fs::File::create(filename).await?;
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
            file.write_all(line.as_bytes()).await?;
        }
        Ok(())
    }

    pub async fn to_pdb(&self, filename: &str) -> std::io::Result<()> {
        let mut file = tokio::fs::File::create(filename).await?;
        let eng = Amber::new(self.sample.system.clone()).energy();
        let pdb = RotateAtDihedral::new(self.sample.system.clone()).to_pdbstring(1, eng);
        file.write_all(pdb.as_bytes()).await?;

        for (i, s) in self.minimized.iter().enumerate() {
            let pdb = RotateAtDihedral::new(s.clone()).to_pdbstring(i + 1, self.energy[i]);
            file.write_all(pdb.as_bytes()).await?;
        }
        Ok(())
    }
}
