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
use clap::ValueEnum;
use nalgebra::{Matrix3, Vector3};
use rand::Rng;
// use rayon::prelude::*;
use std::io::Write;
use std::str::FromStr;
use std::sync::{Arc, LazyLock};

static DIRECTION: LazyLock<Vec<String>> = LazyLock::new(|| {
    include_str!("../library/sobol/new-joe-kuo-6.21201")
        .lines()
        .map(String::from)
        .collect()
});

const SIZE: usize = 32;

pub struct Sobol {
    count: usize,
    total: usize,
    method: Method,
    current: Vec<usize>,
    direction: Vec<Vec<usize>>,
}

impl Sobol {
    pub fn new(dimension: usize, method: Method) -> Self {
        assert!(
            (1..=21201).contains(&dimension),
            "DIMESNSION must in range (1-21201)"
        );
        let mut sobol = Self {
            count: 0,
            total: usize::pow(2, SIZE as u32),
            method,
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

    pub fn get_index(&self, index: usize) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.current.len());

        // For each dimension
        for d in 0..self.current.len() {
            let mut x: usize = 0;
            let mut i: u32 = 0;
            let mut idx = index;

            while idx > 0 {
                if idx & 1 == 1 {
                    x ^= self.direction[d][i as usize];
                }
                idx >>= 1;
                i += 1;
            }

            let numerator = x as f64;
            let denominator = f64::powi(2.0, SIZE as i32);
            result.push(numerator / denominator);
        }
        self.transform_angle(uniform_noise(result.as_slice(), 0.05))
    }

    pub fn transform_angle(&self, angle: Vec<f64>) -> Vec<f64> {
        let scale = 360.0;
        let mut rng = rand::thread_rng();

        match self.method {
            Method::None => angle.iter().map(|x| x * scale).collect(),
            Method::Explore => angle.iter().map(|x| (x * scale) - (scale / 2.0)).collect(),
            Method::Fold => {
                let s = scale / 4.0;
                let transform_idx = rng.gen_range(0..4);
                match transform_idx {
                    0 => angle.iter().map(|x| x * s).collect(),
                    1 => angle.iter().map(|x| (x * s) - s).collect(),
                    2 => angle.iter().map(|x| (x * s) - (2.0 * s)).collect(),
                    _ => angle.iter().map(|x| (x * s) + s).collect(),
                }
            }
            Method::Search => {
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
        }
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
        let mut s = System::from_pdb(file, FF::Amber99SB.init());
        s.particles = s
            .particles
            .iter()
            .filter(|x| x.record == "ATOM")
            .cloned()
            .collect::<Vec<_>>();
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

    pub fn current_dihedral(&self) -> Vec<f64> {
        let mut dihedral = Vec::new();
        for dih in self.system.dihedral.iter() {
            if let Some((a, b, c, d)) = self
                .system
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
        for (a, phi) in self.system.dihedral.iter().zip(angle.clone()) {
            // let phi = theta;
            // let phi = 90.0;
            // if i != 0 || i != self.system.dihedral.len() - 1 {
            //     phi += 180.0
            // }
            if let Some(p) = self.rotated.iter().find(|i| i.serial == a.0.serial) {
                // println!("{:?}", p);
                let v1 = Vector3::new(p.position[0], p.position[1], p.position[2]);
                if let Some(q) = self.rotated.iter().find(|i| i.serial == a.1.serial) {
                    let v2 = Vector3::new(q.position[0], q.position[1], q.position[2]);
                    // println!("V1: {:?}; V2: {:?}, {:?}", v1, v2, angle);
                    let dcos = Self::elemen(v1, v2);
                    for j in self
                        .rotated
                        .iter_mut()
                        .take(a.3.serial)
                        .skip(a.2.serial - 1)
                    {
                        let avector = Vector3::new(j.position[0], j.position[1], j.position[2]);
                        // println!("AVECTOR: {:?} => {}:{}", avector, j.serial, j.name);
                        let r = Self::rotor(dcos, avector - v1, phi) + v1;
                        j.position[0] = r[0];
                        j.position[1] = r[1];
                        j.position[2] = r[2];
                        // println!("NEWPOSITION: {:?}", j.position);
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
        ) * 0.1;
        let u3 = Vector3::new(
            d.position[0] - c.position[0],
            d.position[1] - c.position[1],
            d.position[2] - c.position[2],
        ) * 0.1;
        let sinth = (u2.norm() * u1).dot(&u2.cross(&u3));
        let costh = u1.cross(&u2).dot(&u2.cross(&u3));
        sinth.atan2(costh)
    }

    pub fn energy(&self) -> f64 {
        Amber::new(Arc::clone(&self.system)).energy()
    }

    pub fn to_pdbstring(&self, model: usize, energy: f64) -> String {
        let pdb = parser::atoms_to_pdbstring(self.rotated.clone());
        format!("MODEL{:>9}{:>10.3}\n{}\nENDMDL\n\n", model, energy, pdb)
    }
}

#[derive(Debug, Clone, ValueEnum)]
pub enum Method {
    Fold,
    Search,
    Explore,
    None,
}

impl FromStr for Method {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "fold" => Ok(Method::Fold),
            "search" => Ok(Method::Search),
            "explore" => Ok(Method::Explore),
            "none" => Ok(Method::None),
            _ => Err(format!("Invalid method: {}", s)),
        }
    }
}

#[derive(Clone)]
pub struct Sampler {
    pub method: Method,
    pub folder: String,
    pub energy: Vec<(usize, f64)>,
    pub dihedral: usize,
    pub system: Arc<System>,
    pub rotate: RotateAtDihedral,
}

impl Sampler {
    pub fn new(system: Arc<System>, method: Method, folder: String) -> Self {
        if std::path::Path::new(&folder).exists() {
            std::fs::remove_dir_all(&folder).unwrap();
        }
        std::fs::create_dir_all(folder.to_string()).unwrap();
        Self {
            method,
            folder: folder.to_string(),
            energy: Vec::new(),
            dihedral: system.dihedral.len(),
            system: Arc::clone(&system),
            rotate: RotateAtDihedral::new(Arc::clone(&system)),
        }
    }

    pub fn sample(&mut self, max: usize) {
        // assert!(max > out);

        let s = Sobol::new(self.dihedral, self.method.clone());
        for i in 0..max {
            let mut rotate = RotateAtDihedral::new(Arc::clone(&self.system));
            rotate.rotate(s.get_index(i));
            // let energy_val = rotate.energy();
            // self.energy.push((i, energy_val));

            // Generate the PDB string
            // let pdb_content = rotate.to_pdbstring(i, energy_val);

            // Write to file
            let filename = format!("{}/sample_{:04}.pdb", self.folder, i);
            let mut file = std::fs::File::create(&filename).unwrap();
            file.write_all(rotate.to_pdbstring(i, 0.0).as_bytes())
                .unwrap();

            // Sort self.energy tuples based on energy values (index 1)
            // self.energy.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            // for i in 0..out {
            //     let idx = self.energy[i];

            //     let mut rotate = RotateAtDihedral::new(Arc::clone(&self.system));
            //     // Rotate to that specific configuration
            //     rotate.rotate(s.get_index(idx.0));

            //     // Generate the PDB string
            //     let pdb_content = rotate.to_pdbstring(i, idx.1);

            //     // Write to file
            //     let filename = format!("{}/sample_{:04}.pdb", self.folder, i);
            //     let mut file = std::fs::File::create(&filename).unwrap();
            //     file.write_all(pdb_content.as_bytes()).unwrap();
        }

        let mut file = std::fs::File::create(&format!("{}/sample.out", self.folder)).unwrap();

        let mut energies: Vec<(usize, f64)> = Vec::new();
        for i in 0..max {
            let filename = format!("{}/sample_{:04}.pdb", self.folder, i);
            let mut sys = System::from_pdb(&filename, self.system.forcefield.clone());
            sys.init_parameters();

            let eng = Amber::new(Arc::new(sys)).energy();

            energies.push((i, eng));

            // println!("Energy: {}:  {}kJ/mol", i, eng);
        }

        energies.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        for (i, eng) in energies {
            writeln!(file, "Sample {}: {}kJ/mol", i, eng).unwrap();
        }
    }
}
