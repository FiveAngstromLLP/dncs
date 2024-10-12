// use clap::Parser;
use libdncs::parser::atoms_to_pdbstring;
use libdncs::system::System;
fn main() {
    let s = System::new("YGGFM");
    println!("{}", atoms_to_pdbstring(s.particles))
}
