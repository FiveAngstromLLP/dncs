#![allow(dead_code)]

use crate::forcefield::Amber;
use crate::parser::{self, Atom, FF};
use crate::system::{Particles, System};
use nalgebra::{Matrix3, Vector3};
use rand::Rng;
use rayon::prelude::*;
use std::io::Write;
use std::sync::LazyLock;

static DIRECTION: LazyLock<Vec<String>> = LazyLock::new(|| {
    include_str!("../data/new-joe-kuo-6.21201")
        .lines()
        .map(String::from)
        .collect()
});

const SIZE: usize = 32;

pub struct Sobol {
    count: usize,
    total: usize,
    current: Vec<usize>,
    direction: Vec<Vec<usize>>,
}

impl Sobol {
    pub fn new(dimension: usize) -> Self {
        assert!(
            (1..=21201).contains(&dimension),
            "DIMESNSION must in range (1-21201)"
        );
        let mut sobol = Self {
            count: 0,
            total: usize::pow(2, SIZE as u32),
            current: vec![0; dimension],
            direction: vec![vec![0; SIZE]; dimension],
        };
        sobol.direction = sobol.get_direction(dimension);
        sobol
    }
    fn get_direction(&self, d: usize) -> Vec<Vec<usize>> {
        (1..=d)
            .into_par_iter()
            .map(|d| match d {
                1 => (1..=SIZE).map(|i| 1 << (SIZE - i)).collect(),
                _ => {
                    let mut val = vec![0; SIZE];
                    let direction: Vec<usize> = DIRECTION[d]
                        .split_whitespace()
                        .skip(2)
                        .map(|f| f.parse().unwrap())
                        .collect();
                    for i in 0..SIZE {
                        let s = direction.len() - 1;
                        if i < s {
                            for j in 0..s {
                                val[j] = direction[j + 1] << (SIZE - j - 1)
                            }
                        } else {
                            for k in s..SIZE {
                                val[k] = val[k - s] ^ (val[k - s] >> s);
                                for l in 1..s {
                                    let a = (direction[0] >> (s - l)) & 1;
                                    val[k] ^= a * val[k - l];
                                }
                            }
                        }
                    }
                    val
                }
            })
            .collect()
    }
}

impl Iterator for Sobol {
    type Item = Vec<f64>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.count < self.total {
            let rmz = (self.count ^ (self.total - 1)).trailing_zeros() as usize;
            let val: Vec<usize> = self.current.clone();
            self.count += 1;
            self.current = val
                .par_iter()
                .enumerate()
                .map(|(d, m)| self.direction[d][rmz] ^ m)
                .collect();
            let num: Vec<f64> = val
                .par_iter()
                .map(|i| *i as f64 / f64::powi(2.0, SIZE as i32))
                .collect();
            // let num = bakers_transform(num);
            let num = uniform_noise(&num, 0.05);
            Some(num)
        } else {
            None
        }
    }
}

#[inline]
fn bakers_transform(input: Vec<f64>) -> Vec<f64> {
    let m = input.len();
    if m < 2 {
        input
    } else {
        let k = ((m as f64 * input[0]).floor() as usize + 1).min(m);
        let first = m as f64 * input[0] - (k - 1) as f64;
        let transformed: Vec<f64> = std::iter::once(first)
            .chain(
                input[1..]
                    .iter()
                    .map(|&xi| (xi + (k - 1) as f64) / m as f64),
            )
            .collect();
        transformed
    }
}

#[inline]
fn uniform_noise(point: &[f64], noise_amplitude: f64) -> Vec<f64> {
    let mut rng = rand::thread_rng();
    point
        .iter()
        .map(|&x| {
            let noise = rng.gen_range(-noise_amplitude..=noise_amplitude);
            (x + noise).clamp(0.0, 1.0)
        })
        .collect()
}

#[test]
fn sobol() {
    let s = Sobol::new(10);
    for i in s.take(10) {
        println!("{:?}", i)
    }
}

#[derive(Clone)]
pub struct RotateAtDihedral {
    pub system: System,
    pub rotated: Particles,
}

impl RotateAtDihedral {
    pub fn new(system: System) -> Self {
        Self {
            system: system.clone(),
            rotated: system.particles.clone(),
        }
    }

    pub fn from_pdb(file: &str) -> Vec<f64> {
        let mut s = System::from_pdb(file, FF::AMBER99SB.init());
        let mut scopy = s.clone();
        scopy.get_dihedral();
        s.get_dihedral();

        let mut dihedral = Vec::new();

        for dih in scopy.dihedral {
            if let Some((a, b, c, d)) = s
                .dihedral_angle
                .iter()
                .find(|i| dih.0.serial == i.1.serial && dih.1.serial == i.2.serial)
            {
                // println!(
                //     "{}:{} {}:{} {}:{} {}:{}",
                //     a.name, a.serial, b.name, b.serial, c.name, c.serial, d.serial, d.name,
                // );
                dihedral.push(Self::dihedral_angle(a, b, c, d))
            }
        }

        dihedral

        // s.dihedral
        //     .iter()
        //     .map(|(a, b, c, d)| Self::dihedral_angle(a, b, c, d))
        //     .collect()
    }

    /// Rotate the atoms at Dihedral angle
    pub fn rotate(&mut self, angle: Vec<f64>) {
        for (i, (a, theta)) in self.system.dihedral.iter().zip(angle).enumerate() {
            let mut phi = theta;
            if i != 0 || i != self.system.dihedral.len() - 1 {
                phi += 180.0
            }
            if let Some(p) = self.rotated.iter().find(|i| i.serial == a.0.serial) {
                let v1 = Vector3::new(p.position[0], p.position[1], p.position[2]);
                if let Some(q) = self.rotated.iter().find(|i| i.serial == a.1.serial) {
                    let v2 = Vector3::new(q.position[0], q.position[1], q.position[2]);
                    let dcos = Self::elemen(v1, v2);
                    for j in self
                        .rotated
                        .iter_mut()
                        .take(a.3.serial)
                        .skip(a.2.serial - 1)
                    {
                        let avector = Vector3::new(j.position[0], j.position[1], j.position[2]);
                        let r = Self::rotor(dcos, avector - v1, phi) + v1;
                        j.position[0] = r[0];
                        j.position[1] = r[1];
                        j.position[2] = r[2];
                    }
                }
            }
        }
        self.system.particles = self.rotated.clone();
    }

    pub fn rotated_energy(&mut self, angle: Vec<f64>) -> f64 {
        self.rotate(angle);
        let ff = Amber::new(self.system.clone());
        ff.energy()
    }

    /// Normalization
    #[inline]
    fn elemen(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
        (b - a) / (a - b).norm()
    }

    /// Rotate vector a at the axis of vector b with angle theta
    #[inline]
    fn rotor(a: Vector3<f64>, b: Vector3<f64>, theta: f64) -> Vector3<f64> {
        let phi = theta.to_radians();
        let op1 = Matrix3::<f64>::identity() * phi.cos();
        let op2 = a * a.transpose() * (1.0 - phi.cos());
        let op3 = Matrix3::new(0.0, a[2], -a[1], -a[2], 0.0, a[0], a[1], -a[0], 0.0) * phi.sin();
        let mat = op1 + op2 - op3;
        mat * b
    }

    pub fn dihedral_angle(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> f64 {
        let u1 = Vector3::new(
            b.position[0] - a.position[0],
            b.position[1] - a.position[1],
            b.position[2] - a.position[2],
        );
        let u2 = Vector3::new(
            c.position[0] - b.position[0],
            c.position[1] - b.position[1],
            c.position[2] - b.position[2],
        );
        let u3 = Vector3::new(
            d.position[0] - c.position[0],
            d.position[1] - c.position[1],
            d.position[2] - c.position[2],
        );
        let sinth = (u2.norm() * u1).dot(&u2.cross(&u3));
        let costh = u1.cross(&u2).dot(&u2.cross(&u3));
        sinth.atan2(costh).to_degrees()
    }

    pub fn energy(&self) -> f64 {
        Amber::new(self.system.clone()).energy()
    }

    pub fn to_pdbstring(&self, model: usize, energy: f64) -> String {
        let pdb = parser::atoms_to_pdbstring(self.rotated.clone());
        format!("MODEL{:>9}{:>16.10}\n{}\nENDMDL\n\n", model, energy, pdb)
    }
}

#[derive(Clone)]
pub struct Sampler {
    pub system: System,
    pub rotate: RotateAtDihedral,
    pub angles: Vec<Vec<f64>>,
    pub sample: Vec<System>,
    pub energy: Vec<f64>,
}

impl Sampler {
    pub fn new(system: System) -> Self {
        Self {
            system: system.clone(),
            rotate: RotateAtDihedral::new(system.clone()),
            angles: Vec::new(),
            sample: Vec::new(),
            energy: Vec::new(),
        }
    }

    pub fn sample(&mut self, maxsample: usize) {
        let n = self.system.dihedral.len();
        for phi in Sobol::new(n).skip(32).take(maxsample * 2) {
            let angle: Vec<f64> = phi.iter().map(|i| (i * 360.0) - 180.0).collect();
            self.rotatesample(angle.clone());
            let energy = self.rotate.energy();
            self.energy.push(energy);
            self.angles.push(angle.clone());
            self.sample.push(self.rotate.system.clone());
        }
        self.conformational_sort();
        self.energy = self.energy.clone().into_iter().take(maxsample).collect();
        self.angles = self.angles.clone().into_iter().take(maxsample).collect();
        self.sample = self.sample.clone().into_iter().take(maxsample).collect();
    }

    fn rotatesample(&mut self, angle: Vec<f64>) {
        self.reinit();
        self.rotate.rotate(angle);
    }

    fn conformational_sort(&mut self) {
        const KBT: f64 = 300.0 * 1.380649e-23 * 6.02214076e23 / 4184.0; // KCal/mol
        let weight: Vec<f64> = self.energy.iter().map(|e| (-e / KBT).exp()).collect();
        let z: f64 = weight.iter().sum();
        let normalized: Vec<f64> = weight.iter().map(|w| w / z).collect();
        // println!("Normalized : {normalized:?}");
        let energies = std::mem::take(&mut self.energy);
        let angles = std::mem::take(&mut self.angles);
        let samples = std::mem::take(&mut self.sample);

        let mut combined: Vec<(f64, Vec<f64>, System, f64)> = energies
            .into_iter()
            .zip(angles)
            .zip(samples)
            .zip(normalized)
            .map(|(((e, a), s), w)| (e, a, s, w))
            .collect();

        combined.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));

        self.energy = combined.iter().map(|(e, _, _, _)| *e).collect();
        self.angles = combined.iter().map(|(_, a, _, _)| a.clone()).collect();
        self.sample = combined.into_iter().map(|(_, _, s, _)| s).collect();
    }

    #[inline]
    fn reinit(&mut self) {
        self.rotate = RotateAtDihedral::new(self.system.clone());
    }

    pub fn write_angles(&self, filename: &str) {
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
        let eng = Amber::new(self.system.clone()).energy();
        let pdb = RotateAtDihedral::new(self.system.clone()).to_pdbstring(1, eng);
        file.write_all(pdb.as_bytes()).unwrap();
        for (i, s) in self.sample.iter().enumerate() {
            let pdb = RotateAtDihedral::new(s.clone()).to_pdbstring(i + 1, self.energy[i]);
            file.write_all(pdb.as_bytes()).unwrap();
        }
    }

    pub fn to_pdbfiles(&self, foldername: &str) {
        std::fs::create_dir_all(foldername).unwrap();
        for (i, s) in self.sample.iter().enumerate() {
            let pdb = RotateAtDihedral::new(s.clone()).to_pdbstring(1, self.energy[i]);
            let filename = format!("{}/sample_{:04}.pdb", foldername, i);
            let mut file = std::fs::File::create(filename).unwrap();
            file.write_all(pdb.as_bytes()).unwrap();
        }
    }
}
