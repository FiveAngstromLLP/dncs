#![allow(dead_code)]

use crate::forcefield::Amber;
use crate::parser::{self, Atom};
use crate::system::{Particles, System};
use nalgebra::{Matrix3, Vector3};
use rand::Rng;
use rayon::prelude::*;
use std::sync::{Arc, LazyLock};

const DIRECTION: LazyLock<Vec<String>> = LazyLock::new(|| {
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
            dimension >= 1 && dimension <= 21201,
            "DIMESNSION must in range (1-21201)"
        );
        let mut sobol = Self {
            count: 0,
            total: 1 << SIZE,
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
            let num = bakers_transform(num);
            // let num = uniform_noise(&num, 0.05);
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
        .map(|&x| x + rng.gen_range(-noise_amplitude..=noise_amplitude))
        .collect()
}

#[test]
fn sobol() {
    let s = Sobol::new(10);
    for i in s.take(10) {
        println!("{:?}", i)
    }
}

pub struct RotateAtDihedral {
    pub system: Arc<System>,
    pub rotated: Particles,
}

impl RotateAtDihedral {
    pub fn new(system: Arc<System>) -> Self {
        Self {
            system: Arc::clone(&system),
            rotated: system.particles.clone(),
        }
    }

    /// Rotate the atoms at Dihedral angle
    pub fn rotate(&mut self, angle: Vec<f64>) {
        for (i, (a, theta)) in self.system.dihedral.iter().zip(angle).enumerate() {
            let mut phi = theta;
            if i == 0 || i == self.system.dihedral.len() - 1 {
                phi = phi
            } else {
                phi = phi + 180.0
            }
            if let Some(p) = self.rotated.iter().find(|i| i.serial == a.0) {
                let v1 = Vector3::new(p.position[0], p.position[1], p.position[2]);
                if let Some(q) = self.rotated.iter().find(|i| i.serial == a.1) {
                    let v2 = Vector3::new(q.position[0], q.position[1], q.position[2]);
                    let dcos = Self::elemen(v1, v2);
                    for j in self.rotated.iter_mut().take(a.3).skip(a.2 - 1) {
                        let avector = Vector3::new(j.position[0], j.position[1], j.position[2]);
                        let r = Self::rotor(dcos, avector - v1, phi) + v1;
                        j.position[0] = r[0];
                        j.position[1] = r[1];
                        j.position[2] = r[2];
                    }
                }
            }
        }
    }

    /// Normalization
    #[inline]
    fn elemen(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
        (b - a) / (a - b).norm()
    }

    /// Rotate vector <a> at the axis of vector <b> with angle theta
    #[inline]
    fn rotor(a: Vector3<f64>, b: Vector3<f64>, theta: f64) -> Vector3<f64> {
        let phi = theta.to_radians();
        let op1 = Matrix3::<f64>::identity() * phi.cos();
        let op2 = a * a.transpose() * (1.0 - phi.cos());
        let op3 = Matrix3::new(0.0, a[2], -a[1], -a[2], 0.0, a[0], a[1], -a[0], 0.0) * phi.sin();
        let mat = op1 + op2 - op3;
        mat * b
    }

    fn dihedral_angle(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> f64 {
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
        if u1.cross(&u2).dot(&u2.cross(&u3).cross(&u2)) < 0.0 {
            -sinth.atan2(costh).to_degrees()
        } else {
            sinth.atan2(costh).to_degrees()
        }
    }

    fn to_pdbstring(&self, model: usize) -> String {
        let pdb = parser::atoms_to_pdbstring(self.rotated.clone());
        let energy = Amber::new(Arc::clone(&self.system)).energy();
        format!("MODEL{:>9}{:>16.10}\n{}\nENDMDL\n", model, energy, pdb)
    }
}
