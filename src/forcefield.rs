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
use std::sync::Arc;

use crate::parser::Atom;
use crate::sampling::RotateAtDihedral;
use crate::system::System;
use nalgebra::Vector3;
// use rayon::prelude::*;

pub struct Amber {
    pub system: Arc<System>,
}

impl Amber {
    pub fn new(system: Arc<System>) -> Self {
        Self { system }
    }

    pub fn energy(&self) -> f32 {
        let mut lennard_jones = 0.0;
        let mut electrostatic = 0.0;
        let mut harmonic_bond = 0.0;
        let mut harmonic_angle = 0.0;
        let torsional = self.periodic_torsional_force();

        for iatom in self.system.particles.iter() {
            lennard_jones += self.nb_lennard_jones(iatom);
            electrostatic += self.nb_electrostatic(iatom);
            harmonic_bond += self.harmonic_bond_force(iatom);
            harmonic_angle += self.harmonic_angle_force(iatom);
        }

        // println!("Lennard-Jones energy: {}", lennard_jones);
        // println!("Electrostatic energy: {}", electrostatic);
        // println!("Harmonic bond energy: {}", harmonic_bond);
        // println!("Harmonic angle energy: {}", harmonic_angle);
        // println!("Periodic torsional energy: {}", torsional);

        lennard_jones + electrostatic + harmonic_bond + harmonic_angle + torsional
    }

    #[inline]
    fn nb_electrostatic(&self, iatom: &Atom) -> f32 {
        let mut electrostatic = 0.0;
        for l in self.system.nonbonded[iatom.serial - 1].chunks(2) {
            for jatom in self
                .system
                .particles
                .iter()
                .take(l[1].serial)
                .skip(l[0].serial)
            {
                electrostatic += Self::electrostatic_energy(iatom, jatom);
            }
        }
        for j in self.system.bonded1_4[iatom.serial - 1].iter() {
            if let Some(jatom) = self.system.particles.iter().find(|a| a.serial == j.serial) {
                let nb = self.system.forcefield.nonbonded_force.clone();
                electrostatic +=
                    nb.coulomb14scale as f32 * Self::electrostatic_energy(iatom, jatom);
            }
        }
        // println!(
        //     "Electrostatic energy of {}: {}: {}",
        //     iatom.residue, iatom.name, electrostatic
        // );
        electrostatic // Unit >> kJ/mol
    }

    #[inline]
    fn nb_lennard_jones(&self, iatom: &Atom) -> f32 {
        let mut lennard_jones = 0.0;
        for l in self.system.nonbonded[iatom.serial - 1].chunks(2) {
            for jatom in self
                .system
                .particles
                .iter()
                .take(l[1].serial)
                .skip(l[0].serial)
            {
                lennard_jones += Self::lennard_jones_energy(iatom, jatom);
            }
        }
        for j in self.system.bonded1_4[iatom.serial - 1].iter() {
            if let Some(jatom) = self.system.particles.iter().find(|a| a.serial == j.serial) {
                let nb = self.system.forcefield.nonbonded_force.clone();
                lennard_jones += nb.lj14scale as f32 * Self::lennard_jones_energy(iatom, jatom);
            }
        }
        // println!(
        //     "Lennard-Jones energy of {}: {}: {}",
        //     iatom.residue, iatom.name, lennard_jones
        // );
        lennard_jones // Unit >> kJ/mol
    }

    #[inline]
    fn harmonic_bond_force(&self, iatom: &Atom) -> f32 {
        let hbforce = self.system.forcefield.harmonic_bond_force.bonds.clone();
        let mut energy = 0.0;
        for jatom in self.system.firstbonded[iatom.serial - 1].iter() {
            if let Some(hbf) = hbforce
                .iter()
                .find(|h| h.class1 == iatom.atomtype && h.class2 == jatom.atomtype)
            {
                let d = Self::distance(iatom, jatom);
                let eng = 0.5 * hbf.k as f32 * (d - hbf.length as f32).powi(2);
                energy += eng
            }
        }
        energy // Unit >> kJ/mol
    }

    #[inline]
    fn harmonic_angle_force(&self, iatom: &Atom) -> f32 {
        let haforce = self.system.forcefield.harmonic_angle_force.angles.clone();
        let mut energy = 0.0;
        for jatom in self.system.firstbonded[iatom.serial - 1].iter() {
            for katom in self.system.secondbonded[iatom.serial - 1].iter() {
                if let Some(haf) = haforce.iter().find(|h| {
                    h.class1 == iatom.atomtype
                        && h.class2 == jatom.atomtype
                        && h.class3 == katom.atomtype
                }) {
                    let a = Self::angle(iatom, jatom, katom);

                    energy += 0.5 * haf.k as f32 * (a - haf.angle as f32).powi(2);
                }
            }
        }
        energy // Unit >> kJ/mol
    }

    #[inline]
    fn periodic_torsional_force(&self) -> f32 {
        let periodic = self.system.forcefield.periodic_torsion_force.clone();
        let mut energy = 0.0; // KJ/Mol/nm
        for (a, b, c, d) in self.system.dihedral_angle.iter() {
            let dh = RotateAtDihedral::dihedral_angle(a, b, c, d);

            if let Some(ptf) = periodic.proper.iter().find(|h| {
                h.class1.contains(&a.atomtype)
                    && h.class2.contains(&b.atomtype)
                    && h.class3.contains(&c.atomtype)
                    && h.class4.contains(&d.atomtype)
            }) {
                energy += ptf.k1 as f32
                    * (1.0 + (ptf.periodicity1 as f32 * dh - ptf.phase1 as f32).cos());

                if let (Some(k), Some(n), Some(th)) = (ptf.k2, ptf.periodicity2, ptf.phase2) {
                    energy += k as f32 * (1.0 + (n as f32 * dh - th as f32).cos());
                }

                if let (Some(k), Some(n), Some(th)) = (ptf.k3, ptf.periodicity3, ptf.phase3) {
                    energy += k as f32 * (1.0 + (n as f32 * dh - th as f32).cos());
                }

                if let (Some(k), Some(n), Some(th)) = (ptf.k4, ptf.periodicity4, ptf.phase4) {
                    energy += k as f32 * (1.0 + (n as f32 * dh - th as f32).cos());
                }

                if let (Some(k), Some(n), Some(th)) = (ptf.k5, ptf.periodicity5, ptf.phase5) {
                    energy += k as f32 * (1.0 + (n as f32 * dh - th as f32).cos());
                }

                if let (Some(k), Some(n), Some(th)) = (ptf.k6, ptf.periodicity6, ptf.phase6) {
                    energy += k as f32 * (1.0 + (n as f32 * dh - th as f32).cos());
                }
            }
        }
        energy // Unit >> kJ/mol
    }

    #[inline]
    fn hydrogen_bond_energy(&self) -> f32 {
        self.system
            .hydrogen
            .iter()
            .map(|(iatom, jatom)| {
                let mut hb_energy = 0.0;
                let r = Self::distance(iatom, jatom) * 1e-9;
                hb_energy += 7557.0 * (1.0 / r).powi(12) - 4184.0 * (1.0 / r).powi(10); // Unit >> Kcal/mol
                let hb_ljenergy = Self::lennard_jones_energy(iatom, jatom); // Unit >> (kg.Å^2/s^2)
                let hb_ljenergy = hb_ljenergy * 1e-10_f32.powi(2); // Unit >> (kg.m^2/s^2)
                let hb_ljenergy = hb_ljenergy * 6.02214076e23 / 1000.0; // Unit >> Kcal/mol
                hb_energy - hb_ljenergy
            })
            .sum::<f32>()
    }

    #[inline]
    fn lennard_jones_energy(i: &Atom, j: &Atom) -> f32 {
        let r = Self::distance(i, j); // Unit >> nm
        if r == 0.0 {
            return 0.0; // Avoid division by zero
        }
        let sigma = ((i.sigma + j.sigma) / 2.0) as f32; // Unit >> nm
        let epsilon = (i.epsilon * j.epsilon).sqrt() as f32; // Unit >> kJ/mol
        4.0 * epsilon as f32 * ((sigma / r).powi(12) - (sigma / r).powi(6)) // Unit >> kJ/mol
    }

    #[inline]
    fn electrostatic_energy(i: &Atom, j: &Atom) -> f32 {
        let r = Self::distance(i, j); // Unit >> nm
        let q1 = i.charge as f32 * 1.602176634e-19;
        let q2 = j.charge as f32 * 1.602176634e-19;
        let k = (8.9875517682e9 * 6.02214076e23) / 1000.0; // Coulomb's constant in vacuum (kJ/mol)
        k * q1 * q2 / r // Unit >> kJ/mol
    }

    #[inline]
    fn angle(i: &Atom, j: &Atom, k: &Atom) -> f32 {
        let rij = Vector3::new(
            j.position[0] - i.position[0],
            j.position[1] - i.position[1],
            j.position[2] - i.position[2],
        ) * 0.1; // Unit >> nm
        let rkj = Vector3::new(
            j.position[0] - k.position[0],
            j.position[1] - k.position[1],
            j.position[2] - k.position[2],
        ) * 0.1; // Unit >> nm

        let cos_theta = rij.dot(&rkj) / (rij.magnitude() * rkj.magnitude());
        cos_theta.acos() as f32
    }

    #[inline]
    fn distance(i: &Atom, j: &Atom) -> f32 {
        let rij = Vector3::new(
            j.position[0] - i.position[0],
            j.position[1] - i.position[1],
            j.position[2] - i.position[2],
        );
        rij.map(|r| r.powi(2)).map(|x| x as f32).sum().sqrt() * 0.1 // Unit >> nm
    }
}
