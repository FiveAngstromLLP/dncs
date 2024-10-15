extern crate libdncs;

use libdncs::{sampling::RotateAtDihedral, system::System};

fn main() {
    let mut s = System::new("YGGFM");
    s.get_dihedralatoms(true);
    println!("{:?}", s.dihedral.len());
    let y = RotateAtDihedral::from_pdb("examples/test.pdb", true);
    println!("{:?}", y.len());
}
