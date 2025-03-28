use std::sync::Arc;

use libdncs::*;

// Configuration
const FOLDER: &str = "Result/DNCS";
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AmberFB15;
const NO_OF_SAMPLE: usize = 100;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // System
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    // Sample
    let mut sample = Sampler::new(Arc::new(sys), Method::Fold, FOLDER.to_string());
    println!("Generating Samples..");
    sample.sample(NO_OF_SAMPLE, 300.0);
    // Minimizer
    // let mut minimizer = Minimizer::new(FOLDER.to_string(), FORCE_FIELD.init(), 10);
    println!("Executing Minimizer..");
    // minimizer.minimize_all();
    println!("Completed!");
    Ok(())
}
