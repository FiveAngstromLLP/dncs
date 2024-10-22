extern crate libdncs;

use libdncs::forcefield::Amber;
use libdncs::parser::AMBER99SB;
use libdncs::system::System;

fn main() {
    let mut polymer = System::new("YGGFM", (*AMBER99SB).clone());
    polymer.init_parameters();
    polymer.get_dihedral();
    for (a, b, c, d) in polymer.dihedral.iter() {
        println!("{} {} {} {}", a.name, b.name, c.name, d.name)
    }
    let ff = Amber::new(polymer);
    println!("Energy: {}", ff.energy());
}
