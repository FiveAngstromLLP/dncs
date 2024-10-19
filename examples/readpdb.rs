extern crate libdncs;

use libdncs::parser::atoms_to_pdbstring;
use libdncs::sampling::{RotateAtDihedral, Sampler};
use libdncs::system::System;

const S: bool = true;

fn main() {
    let mut system = System::new("YG"); // Sequence
    system.init_parameters();
    system.get_dihedralatoms(S); // include_sidechain
    let mut sample = Sampler::new(system);
    for (e, i) in sample.system.dihedral.iter().enumerate() {
        println!(
            "{} {}:{} {}:{} {}:{} {}:{}",
            e + 1,
            i.0.name,
            i.0.serial,
            i.1.name,
            i.1.serial,
            i.2.name,
            i.2.serial,
            i.3.name,
            i.3.serial
        )
    }
    println!("\n");
    sample.sample(1);
    println!("{:?}", sample.sample.len());

    let string = atoms_to_pdbstring(sample.sample[0].particles.clone());
    std::fs::write("test.pdb", string).unwrap();
    let string = atoms_to_pdbstring(sample.system.particles);
    std::fs::write("system.pdb", string).unwrap();

    let angle = RotateAtDihedral::from_pdb("test.pdb", S);

    println!("{:?}", sample.angles[0].clone());
    println!("{:6.3?}", angle);
}
