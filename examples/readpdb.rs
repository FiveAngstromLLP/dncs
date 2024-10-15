extern crate libdncs;

use libdncs::sampling::Sampler;
use libdncs::system::System;

fn main() {
    let mut system = System::new("YGGFM");
    system.get_dihedralatoms(true);
    let mut sample = Sampler::new(system);
    sample.sample(10);
    println!("{:?}", sample.angles)
}
