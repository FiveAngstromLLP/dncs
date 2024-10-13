mod parser;
mod system;

use crate::parser::atoms_to_pdbstring;
use crate::system::System;
fn main() {
    let s = System::new("YGGFM");
    println!("{}", atoms_to_pdbstring(s.particles))
}
