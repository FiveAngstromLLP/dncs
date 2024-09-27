// use clap::Parser;
use libdncs::parser::generate;
fn main() {
    let s = generate("YGGFM");
    for atom in s {
        println!("{:?}", atom)
    }
}
