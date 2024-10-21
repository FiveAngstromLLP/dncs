extern crate libdncs;

use libdncs::parser::AMBER99SB;
use libdncs::system::System;

fn main() {
    let sys = System::new("ACDEFGHIKLMNPQRSTVWY", (*AMBER99SB).clone());
    for i in &sys.particles {
        println!("{}: {:?}; {:?}", i.residue, i.name, sys.get_atomtype(i));
    }
}
