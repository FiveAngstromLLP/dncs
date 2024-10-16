extern crate libdncs;

use libdncs::sampling::Sampler;
use libdncs::system::System;

fn main() {
    let mut system = System::new("YGGFM"); // Sequence
    system.init_parameters();
    system.get_dihedralatoms(true); // include_sidechain
    let mut sample = Sampler::new(system);
    sample.sample(10); // Samples
    sample.conformational_sort();
    sample.write_sampled_angles("angles.out");
    sample.to_pdb("Sample.pdb")
}
