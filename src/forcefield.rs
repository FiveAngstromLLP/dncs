#![allow(dead_code)]
use crate::parser::Atom;
use crate::system::System;
use nalgebra::Vector3;
use rayon::prelude::*;

pub struct Amber {
    pub system: System,
}

impl Amber {
    pub fn new(system: System) -> Self {
        Self { system }
    }

    pub fn energy(&self) -> f64 {
        let energy = self
            .system
            .particles
            .par_iter()
            .map(|iatom| self.nonbonded_energy(iatom)) // Unit >> (kg.Å^2/s^2)
            .map(|energy| energy * 1e-10_f64.powi(2)) // Unit >> (kg.m^2/s^2)
            .map(|energy| energy * 6.02214076e23 / 4184.0) // Unit >> Kcal/mol
            .sum::<f64>();
        energy + self.hydrogen_bond_energy()
    }

    pub fn potential(&mut self) {
        let potentials: Vec<_> = self
            .system
            .particles
            .par_iter()
            .map(|iatom| self.nonbonded_potential(iatom)) // Unit >> (kg.Å/s^2)
            .collect();
        for (iatom, potential) in self.system.particles.iter_mut().zip(potentials) {
            iatom.force[0] = potential.x; // Unit >> (kg.Å/s^2)
            iatom.force[1] = potential.y; // Unit >> (kg.Å/s^2)
            iatom.force[2] = potential.z; // Unit >> (kg.Å/s^2)
        }
    }

    fn nonbonded_energy(&self, iatom: &Atom) -> f64 {
        let mut lennard_jones = 0.0;
        let mut electrostatic = 0.0;
        for l in self.system.nonbonded[iatom.serial - 1].chunks(2) {
            for jatom in self
                .system
                .particles
                .iter()
                .take(l[1].serial)
                .skip(l[0].serial)
            {
                lennard_jones += Self::lennard_jones_energy(iatom, jatom);
                electrostatic += Self::electrostatic_energy(iatom, jatom);
            }
        }
        for j in self.system.bonded1_4[iatom.serial - 1].iter() {
            if let Some(jatom) = self.system.particles.iter().find(|a| a.serial == j.serial) {
                lennard_jones += 0.500 * Self::lennard_jones_energy(iatom, jatom);
                electrostatic += 0.500 * Self::electrostatic_energy(iatom, jatom);
            }
        }
        lennard_jones + electrostatic // Unit >> (kg.Å^2/s^2)
    }

    fn nonbonded_potential(&self, iatom: &Atom) -> Vector3<f64> {
        let mut lennard_jones = Vector3::zeros();
        let mut electrostatic = Vector3::zeros();
        for l in self.system.nonbonded[iatom.serial - 1].chunks(2) {
            for jatom in self
                .system
                .particles
                .iter()
                .take(l[1].serial)
                .skip(l[0].serial)
            {
                lennard_jones += Self::lennard_jones_potential(iatom, jatom);
                electrostatic += Self::electrostatic_potential(iatom, jatom);
            }
        }
        for j in self.system.bonded1_4[iatom.serial - 1].iter() {
            if let Some(jatom) = self.system.particles.iter().find(|a| a.serial == j.serial) {
                lennard_jones += 0.500 * Self::lennard_jones_potential(iatom, jatom);
                electrostatic += 0.500 * Self::electrostatic_potential(iatom, jatom);
            }
        }
        lennard_jones + electrostatic // Unit >> (kg.Å/s^2)
    }

    // fn harmonic_bond_force(&self, iatom: &Atom) -> f64 {

    //     if let Some(atype) = self.system.get_atomtype(iatom) {
    //         let id = atype.;

    //     }
    //     0.0
    // }

    fn hydrogen_bond_energy(&self) -> f64 {
        self.system
            .hydrogen
            .iter()
            .map(|(iatom, jatom)| {
                let mut hb_energy = 0.0;
                let r = Self::distance(iatom, jatom);
                hb_energy += 7557.0 * (1.0 / r).powi(12) - 4184.0 * (1.0 / r).powi(10); // Unit >> Kcal/mol
                let hb_ljenergy = Self::lennard_jones_energy(iatom, jatom); // Unit >> (kg.Å^2/s^2)
                let hb_ljenergy = hb_ljenergy * 1e-10_f64.powi(2); // Unit >> (kg.m^2/s^2)
                let hb_ljenergy = hb_ljenergy * 6.02214076e23 / 4184.0; // // Unit >> Kcal/mol
                hb_energy - hb_ljenergy
            })
            .sum::<f64>()
    }

    #[inline]
    fn lennard_jones_energy(i: &Atom, j: &Atom) -> f64 {
        let r = Self::distance(i, j); // Unit >> Å
        let sigma = (i.sigma + j.sigma) / 2.0; // Unit >> Å
        let epsilon = (i.epsilon * j.epsilon).sqrt(); // Unit >> kcal/mol
        let epsilon = epsilon * 4184.0 / 6.02214076e23; // Unit >> Joules (or) Kg.m^2/s^2
        let epsilon = epsilon * 1e10_f64.powi(2); // Unit >> (kg.Å^2/s^2)
        4.0 * epsilon * ((sigma / r).powi(12) - (sigma / r).powi(6)) // Unit >> (kg.Å^2/s^2)
    }

    #[inline]
    fn lennard_jones_potential(i: &Atom, j: &Atom) -> Vector3<f64> {
        let r = Self::distance(i, j); // Unit >> Å
        let sigma = (i.sigma + j.sigma) / 2.0; // Unit >> Å
        let epsilon = (i.epsilon * j.epsilon).sqrt(); // Unit >> kcal/mol
        let epsilon = epsilon * 4184.0 / 6.02214076e23; // Unit >> Joules (or) Kg.m^2/s^2
        let force = -4.0
            * epsilon * 1e10_f64.powi(2) // Epsilon Unit >> (kg.Å^2/s^2)
            * ((12.0 * sigma.powi(12) / r.powi(13)) - (6.0 * sigma.powi(6) / r.powi(7))); // Unit >> (1/Å)
        force * Self::vector_distance(i, j) / r // Unit >> (kg.Å/s^2)
    }

    #[inline]
    fn electrostatic_energy(i: &Atom, j: &Atom) -> f64 {
        let r = Self::distance(i, j); // Unit >> Å
        let q1 = i.charge * 1.602176634e-19f64; // Unit >> C
        let q2 = j.charge * 1.602176634e-19f64; // Unit >> C
        let k = 8.9875517923e9 * 1e10_f64.powi(3); // Unit >> (Kg.Å^3)/(s^2.C^2)
        k * q1 * q2 / r // Unit >> (kg.Å^2/s^2)
    }

    #[inline]
    fn electrostatic_potential(i: &Atom, j: &Atom) -> Vector3<f64> {
        let r = Self::distance(i, j); // Unit >> Å
        let q1 = i.charge * 1.602176634e-19f64; // Unit >> C
        let q2 = j.charge * 1.602176634e-19f64; // Unit >> C
        let k = 8.9875517923e9 * 1e10_f64.powi(3); // Unit >> (Kg.Å^3)/(s^2.C^2)
        let force = -k * q1 * q2 / r.powi(2); // Unit >> Kg.Å/s^2
        force * Self::vector_distance(i, j) / r // Unit >> Kg.Å/s^2
    }

    #[inline]
    fn distance(i: &Atom, j: &Atom) -> f64 {
        let rij = Vector3::new(
            j.position[0] - i.position[0],
            j.position[1] - i.position[1],
            j.position[2] - i.position[2],
        );
        rij.map(|r| r.powi(2)).sum().sqrt() // Unit >> Å
    }

    #[inline]
    fn vector_distance(i: &Atom, j: &Atom) -> Vector3<f64> {
        Vector3::new(
            j.position[0] - i.position[0],
            j.position[1] - i.position[1],
            j.position[2] - i.position[2],
        ) // Unit >> Å
    }
}
