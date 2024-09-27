// use clap::Parser;
use libdncs::system::System;
fn main() {
    let s = System::new("YGGFM");
    for atom in s.particles {
        println!("{:?}", atom)
    }
}
