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
#![allow(dead_code)]

use crate::forcefield::Amber;
use crate::parser::{self, Atom, FF};
use crate::system::{Particles, System};
use nalgebra::{Matrix3, Vector3};
use rand::Rng;
// use rayon::prelude::*;
use std::io::Write;
use std::sync::{Arc, LazyLock};

static DIRECTION: LazyLock<Vec<String>> = LazyLock::new(|| {
    include_str!("../data/sobol/new-joe-kuo-6.21201")
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
            .into_iter()
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
                .iter()
                .enumerate()
                .map(|(d, m)| self.direction[d][rmz] ^ m)
                .collect();
            let num: Vec<f64> = val
                .iter()
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
    pub system: Arc<System>,
    pub rotated: Particles,
}

impl RotateAtDihedral {
    pub fn new(system: Arc<System>) -> Self {
        Self {
            system: Arc::clone(&system),
            rotated: system.particles.to_vec(),
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
                dihedral.push(Self::dihedral_angle(a, b, c, d))
            }
        }
        dihedral
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
        let sys = Arc::make_mut(&mut self.system);
        sys.particles = self.rotated.clone();
    }

    pub fn rotated_energy(&mut self, angle: Vec<f64>) -> f64 {
        self.rotate(angle);
        let ff = Amber::new(Arc::clone(&self.system));
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
        Amber::new(Arc::clone(&self.system)).energy()
    }

    pub fn to_pdbstring(&self, model: usize, energy: f64) -> String {
        let pdb = parser::atoms_to_pdbstring(self.rotated.clone());
        format!("MODEL{:>9}{:>16.10}\n{}\nENDMDL\n\n", model, energy, pdb)
    }
}

#[derive(Clone)]
pub struct Sampler {
    pub system: Arc<System>,
    pub rotate: RotateAtDihedral,
    pub grid: usize,
    pub angles: Vec<Vec<f64>>,
    pub sample: Vec<Arc<System>>,
    pub energy: Vec<f64>,
}

impl Sampler {
    pub fn new(system: Arc<System>, grid: usize) -> Self {
        Self {
            system: Arc::clone(&system),
            rotate: RotateAtDihedral::new(Arc::clone(&system)),
            grid,
            angles: Vec::new(),
            sample: Vec::new(),
            energy: Vec::new(),
        }
    }

    pub fn sample(&mut self, maxsample: usize) {
        let n = self.system.dihedral.len();
        for phi in Sobol::new(n).skip(32).take(2048) {
            let angle: Vec<f64> = self.transform_angle(phi);
            self.rotatesample(angle.clone());
            let energy = self.rotate.energy();
            self.energy.push(energy);
            self.angles.push(angle.clone());
            self.sample.push(Arc::clone(&self.rotate.system));
            println!("Sampling {}/2048", self.angles.len());
        }
        self.conformational_sort();
        self.energy = self.energy.iter().take(maxsample).map(|e| *e).collect();
        self.angles = self
            .angles
            .iter()
            .take(maxsample)
            .map(|a| a.clone())
            .collect();
        self.sample = self
            .sample
            .iter()
            .take(maxsample)
            .map(|s| Arc::clone(&s))
            .collect();
    }

    pub fn transform_angle(&self, angle: Vec<f64>) -> Vec<f64> {
        let scale = 360.0;
        let mut rng = rand::thread_rng();

        match self.grid {
            1 => angle.iter().map(|x| (x * scale) - (scale / 2.0)).collect(),
            2 => {
                let s = scale / self.grid as f64;
                let transform_idx = rng.gen_range(0..2);
                match transform_idx {
                    0 => angle.iter().map(|x| x * s).collect(),
                    _ => angle.iter().map(|x| (x * s) - s).collect(),
                }
            }
            4 => {
                let s = scale / self.grid as f64;
                let transform_idx = rng.gen_range(0..4);
                match transform_idx {
                    0 => angle.iter().map(|x| x * s).collect(),
                    1 => angle.iter().map(|x| (x * s) - s).collect(),
                    2 => angle.iter().map(|x| (x * s) - (2.0 * s)).collect(),
                    _ => angle.iter().map(|x| (x * s) + s).collect(),
                }
            }
            5 => {
                let s1 = scale / 2.0;
                let s = scale / 4.0;
                let transform_idx = rng.gen_range(0..7);
                match transform_idx {
                    0 => angle.iter().map(|x| (x * scale) - (scale / 2.0)).collect(),
                    1 => angle.iter().map(|x| x * s1).collect(),
                    2 => angle.iter().map(|x| (x * s1) - s1).collect(),
                    3 => angle.iter().map(|x| x * s).collect(),
                    4 => angle.iter().map(|x| (x * s) - s).collect(),
                    5 => angle.iter().map(|x| (x * s) - (2.0 * s)).collect(),
                    _ => angle.iter().map(|x| (x * s) + s).collect(),
                }
            }
            _ => angle,
        }
    }

    fn rotatesample(&mut self, angle: Vec<f64>) {
        self.reinit();
        self.rotate.rotate(angle);
    }

    fn conformational_sort(&mut self) {
        // const KBT: f64 = 300.0 * 1.380649e-23 * 6.02214076e23 / 4184.0; // KCal/mol

        // let weight: Vec<f64> = self.energy.iter().map(|e| (-e / KBT).exp()).collect();
        // let z: f64 = weight.iter().sum();
        // let normalized: Vec<f64> = weight.iter().map(|w| w / z).collect();

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
        let sorted_sample: Vec<Arc<System>> = indices
            .iter()
            .map(|&i| Arc::clone(&self.sample[i]))
            .collect();

        // Replace original vectors with sorted ones
        let minenergy = sorted_energy.iter().copied().fold(f64::INFINITY, f64::min);
        self.energy = sorted_energy.iter().map(|eng| eng - minenergy).collect();
        self.angles = sorted_angles;
        self.sample = sorted_sample;
    }

    #[inline]
    fn reinit(&mut self) {
        self.rotate = RotateAtDihedral::new(Arc::clone(&self.system));
    }

    pub fn write_angles(&self, filename: &str) {
        let mut file = std::fs::File::create(filename).unwrap();
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
            file.write_all(line.as_bytes()).unwrap();
        }
    }

    pub fn to_pdb(&self, filename: &str) {
        let mut file = std::fs::File::create(filename).unwrap();

        for (i, s) in self.sample.iter().enumerate() {
            let pdb = RotateAtDihedral::new(Arc::clone(s)).to_pdbstring(i + 1, self.energy[i]);
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

        let mut file = std::fs::File::create(format!("{}/../structure.pdb", foldername)).unwrap();
        let eng = Amber::new(Arc::clone(&self.system)).energy();
        let pdb = RotateAtDihedral::new(Arc::clone(&self.system)).to_pdbstring(1, eng);
        file.write_all(pdb.as_bytes()).unwrap();
    }
}
