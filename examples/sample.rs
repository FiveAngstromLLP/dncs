use std::sync::Arc;

use libdncs::*;

// Configuration
const NAME: &str = "DNCS";
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AmberFB15;
const NO_OF_SAMPLE: usize = 10;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let folder = format!("Result/{}", NAME);
    let mut sample = Sampler::new(Arc::new(sys), Method::Fold, folder.to_string());
    sample.sample(NO_OF_SAMPLE, 300.0);
}
