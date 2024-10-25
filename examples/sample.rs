use libdncs::*;

// Configuration
const NAME: &str = "DNCS";
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;
const NO_OF_SAMPLE: usize = 10;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let mut sample = Sampler::new(sys);
    sample.sample(NO_OF_SAMPLE);
    sample.write_angles(&format!("{}.out", NAME));
    sample.to_pdb(&format!("{}.pdb", NAME));
}
