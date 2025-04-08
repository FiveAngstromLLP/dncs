use std::sync::Arc;

use libdncs::*;

// Configuration

const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AmberFB15;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.to_pdb("sample.pdb");
    sys.init_parameters();

    let mut sample = Sampler::new(Arc::new(sys), Method::Search, "Result".to_string());
    sample.sample(10000);
}
