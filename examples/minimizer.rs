use libdncs::*;

// Configuration
const NAME: &str = "DNCS";
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;
const NO_OF_SAMPLE: usize = 10;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // System
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    // Sample
    let mut sample = Sampler::new(sys);
    println!("Generating Samples..");
    sample.sample(NO_OF_SAMPLE);
    // Minimizer
    let mut minimizer = Minimizer::new(sample);
    println!("Executing Minimizer..");
    minimizer.minimize();
    minimizer.conformational_sort();
    minimizer.write_angles(&format!("{}.out", NAME))?;
    minimizer.to_pdb("minimized.pdb")?;
    println!("Completed!");
    Ok(())
}
