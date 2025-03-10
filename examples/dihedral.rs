extern crate libdncs;

use libdncs::{System, FF};

fn main() {
    let mut sys = System::new("YGGFM", FF::AmberFB15.init());
    sys.init_parameters();
    for d in sys.dihedral {
        println!(
            "{} {} {} {}",
            d.0.serial, d.1.serial, d.2.serial, d.3.serial
        );
    }
}
