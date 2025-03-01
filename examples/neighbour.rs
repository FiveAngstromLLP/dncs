extern crate libdncs;

use libdncs::{System, FF};

fn main() {
    let mut system = System::new("YGGFM", FF::Amber99SB.init());
    system.get_neighbours();

    for (index, neighbour) in system.nonbonded.iter().enumerate() {
        print!(" {}::>>", index + 1);
        for i in neighbour {
            print!("\t{}:{}\t", i.name, i.serial);
        }
        println!();
    }
}
