extern crate libdncs;

use libdncs::{parser::FF, system::System};

fn main() {
    let sys = System::new("ACDEFGHIKLMNPQRSTVWY", FF::AMBER99SB.init());
    for i in &sys.particles {
        println!("{}: {:?}", i.residue, i.name);
    }
}
