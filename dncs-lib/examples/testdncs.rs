extern crate libdncs;

// use libdncs::forcefield::Amber;
use libdncs::sampling::Sampler;
use libdncs::system::System;

fn main() {
    let mut polymer = System::new("YGGFM");
    polymer.init_parameters();
    let mut sample = Sampler::new(polymer.clone());
    sample.system.get_dihedralatoms(true);
    sample.sample(10);
    sample.to_pdb("sample.pdb");
}
