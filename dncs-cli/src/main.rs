// use clap::Parser;
use libdncs::parser::Polygen;
fn main() {
    let s = Polygen::new("YGGFM");
    for atom in s.atoms {
        println!("{:?}", atom)
    }
}
