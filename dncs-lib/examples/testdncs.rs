extern crate libdncs;

// use libdncs::forcefield::Amber;
use libdncs::sampling::Sampler;
use libdncs::system::System;

fn main() {
    let mut polymer = System::new("YGGFM");
    polymer.init_parameters();
    let mut sample = Sampler::new(polymer.clone());
    sample.system.get_dihedralatoms(false);
    sample.sample(100);
    sample.conformational_sort();
    sample.to_pdb("sample.pdb");
}
