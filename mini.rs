#[derive(Clone)]
pub struct Minimizer {
    data: InputParameters,
    pub sample: Sampling,
    pub minized: Sampling,
}

impl Minimizer {
    /// New Minimizer
    pub fn new(data: InputParameters, intake: usize) -> Self {
        let mut s = Sampling::new(data.clone());
        s.generate(intake);
        s.filter_weight();
        Self {
            data: data.clone(),
            sample: s.clone(),
            minized: Sampling::new(data),
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
        for (angle, _) in self.sample.sample.clone() {
            let evaluate = |x: &[f64], gx: &mut [f64]| {
                let val = Self::diff(
                    Box::new(|x: Vec<f64>| self.sample.clone().calculate_for_angle(x)),
                    x.to_vec(),
                );
                gx.copy_from_slice(&val.1);
                Ok(val.0)
            };
            let mut theta: Vec<f64> = angle.clone();
            let prbs = lbfgs()
                .with_max_iterations(10)
                .minimize(&mut theta, evaluate, |prgr| {
                    println!("iter: {:?}", prgr.niter);
                    println!("x: {:?}", prgr.x);
                    println!("Energy: {:?}", prgr.fx);
                    false
                });
            match prbs {
                Ok(p) => {
                    println!("Minimized theta = {:?}", theta);
                    self.minized.sample_for_angle(theta);
                    println!("Minimized Energy = {:?}", p.fx);
                }
                Err(e) => {
                    println!("{:?}", e);
                }
            }
        }
    }

    /// Write Sample Energy
    pub fn write_sample(&mut self) {
        self.sample.write(&format!("{}_sample", self.data.filename))
    }

    /// Write Minimized Energy
    pub fn write_mini(&mut self) {
        self.minized.write(&format!("{}_mini", self.data.filename))
    }
}

// -----------------------------------------------------------------------------
//
use lbfgs::*;
use rayon::prelude::*;
use std::sync::Arc;

pub struct Minimizer {
    sample: Arc<Sampler>,
    minimized: Vec<(f64, System)>,
}

impl Minimizer {
    pub fn new(sample: Sampler) -> Self {
        Self {
            sample: Arc::new(sample),
            minimized: Vec::new(),
        }
    }

    pub fn minimize(&mut self) {
        let minimized: Vec<(f64, System)> = self
            .sample
            .par_iter()
            .map(|system| {
                let mut lbfgs = Self::create_lbfgs(&system);
                Self::minimize_system(system, &mut lbfgs)
            })
            .collect();

        self.minimized = minimized;
    }

    fn create_lbfgs(system: &System) -> Lbfgs {
        let problem_size = 3 * system.num_atoms();
        let lbfgs_memory_size = 5;
        Lbfgs::new(problem_size, lbfgs_memory_size)
            .with_sy_epsilon(1e-8)
            .with_cbfgs_alpha(1.0)
            .with_cbfgs_epsilon(1e-4)
    }

    fn minimize_system(mut system: System, lbfgs: &mut Lbfgs) -> (f64, System) {
        let mut x = system.get_coordinates();
        let mut g = system.get_gradient();
        let mut energy = system.get_energy();

        let max_iterations = 100;
        for _ in 0..max_iterations {
            lbfgs.apply_hessian(&mut g);

            // Simplified line search
            let step_size = 0.01;
            let x_new: Vec<f64> = x
                .iter()
                .zip(g.iter())
                .map(|(xi, gi)| xi - step_size * gi)
                .collect();

            system.set_coordinates(&x_new);
            let new_energy = system.get_energy();
            let new_g = system.get_gradient();

            let s: Vec<f64> = x_new
                .iter()
                .zip(x.iter())
                .map(|(new, old)| new - old)
                .collect();
            let y: Vec<f64> = new_g
                .iter()
                .zip(g.iter())
                .map(|(new, old)| new - old)
                .collect();

            if lbfgs.update_hessian(&s, &y) == UpdateStatus::UpdateOk {
                x = x_new;
                g = new_g;
                energy = new_energy;
            } else {
                break;
            }

            if g.iter().all(|&gi| gi.abs() < 1e-6) {
                break;
            }
        }

        (energy, system)
    }

    pub fn to_pdb(&self, filename: &str) {
        let mut file = std::fs::File::create(filename).unwrap();
        for (i, (energy, sys)) in self.minimized.iter().enumerate() {
            let pdb = RotateAtDihedral::new(sys.clone()).to_pdbstring(i + 1, *energy);
            file.write_all(pdb.as_bytes()).unwrap();
        }
    }

    pub fn to_pdbfiles(&self, foldername: &str) {
        std::fs::create_dir_all(foldername).unwrap();
        for (i, (energy, sys)) in self.minimized.iter().enumerate() {
            let pdb = RotateAtDihedral::new(sys.clone()).to_pdbstring(i + 1, *energy);
            let filename = format!("{}/minimized_{:04}.pdb", foldername, i);
            let mut file = std::fs::File::create(filename).unwrap();
            file.write_all(pdb.as_bytes()).unwrap();
        }
    }
}
