mod forcefield;
mod minimizer;
mod parser;
mod sampling;
mod system;

use crate::minimizer::Minimizer;
use crate::sampling::Sampler;
use crate::system::System;
use clap::{Arg, Command};
use serde_json;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("DNCS")
        .version("0.1.3")
        .about("Digital Nets Conformational Sample (DNCS)")
        .arg(
            Arg::new("config")
                .short('c')
                .long("config")
                .help("Generates dncs.json file")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    if matches.get_one::<bool>("config").is_some_and(|&v| v) {
        generate_config_file()?;
        return Ok(());
    }

    // Read config from JSON file
    let config_content = match fs::read_to_string("dncs.json") {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Error reading dncs.json: {}", e);
            std::process::exit(1);
        }
    };
    let config: serde_json::Value = match serde_json::from_str(&config_content) {
        Ok(value) => value,
        Err(e) => {
            eprintln!("Failed to parse JSON: {}", e);
            std::process::exit(1);
        }
    };

    let generate = match config["Generate"].as_object() {
        Some(object) => object,
        None => {
            eprintln!("Failed to get Generate object from config");
            std::process::exit(1);
        }
    };
    let molecule = generate["molecule"].as_str().unwrap_or("Sample");
    let sequence = generate["sequence"].as_str().unwrap_or("YGGFM");
    let include_sidechain = generate["include_sidechain"].as_bool().unwrap_or(false);
    let samples = generate["n_samples"].as_u64().unwrap_or(10) as usize;

    let mut system = System::new(sequence);
    system.init_parameters();
    system.get_dihedralatoms(include_sidechain);
    let mut sample = Sampler::new(system);
    sample.sample(samples);
    sample.conformational_sort();

    let result_dir = format!("Result/{}", molecule);
    if fs::metadata(&result_dir).is_ok() {
        fs::remove_dir_all(&result_dir)?;
    }
    fs::create_dir_all(&result_dir)?;

    // sample.write_sampled_angles(&format!("Result/{}/Sampled.out", molecule));
    // sample.to_pdb(&format!("Result/{}/Sampled.pdb", molecule));

    let mut mini = Minimizer::new(sample);
    mini.minimize();
    mini.conformational_sort();
    mini.write_sampled_angles(&format!("Result/{}/{}.out", molecule, molecule));
    mini.to_pdb(&format!("Result/{}/{}.pdb", molecule, molecule));
    Ok(())
}

fn generate_config_file() -> Result<(), std::io::Error> {
    let config_content = r#"
{
    "Generate": {
        "molecule": "Sample",
        "sequence": "AAAAA",
        "include_sidechain": true,
        "n_samples": 10
    }
}
"#;
    fs::write("dncs.json", config_content)?;
    println!("Generated dncs.json file");
    Ok(())
}
