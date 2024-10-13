extern crate libdncs;

use libdncs::forcefield::Amber;
use libdncs::system::System;

fn main() {
    let mut polymer = System::new("YGGFM");
    polymer.init_parameters();
    polymer.to_pdb("YGGFM.pdb");
    let ff = Amber::new(polymer);
    println!("Energy: {}", ff.energy());
}
