extern crate libdncs;

use std::sync::Arc;

use libdncs::forcefield::Amber;
use libdncs::system::System;

fn main() {
    let mut polymer = System::new("YGGFM");
    polymer.init_parameters();
    let amber = Amber::new(Arc::new(polymer));
    println!("Energy: {}", amber.energy());
}
