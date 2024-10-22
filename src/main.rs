mod forcefield;
mod minimizer;
mod parser;
mod sampling;
mod system;
use crate::minimizer::Minimizer;
use crate::sampling::Sampler;
use crate::system::System;
use clap::{Arg, Command};
use parser::FF;
use serde_json::{json, Value};
use std::fs;
use std::path::Path;

struct SimulationParams {
    molecule: String,
    sequence: String,
    n_samples: usize,
    forcefield: String,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
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
        .arg(
            Arg::new("molecule")
                .short('m')
                .long("molecule")
                .help("Molecule name")
                .value_name("NAME"),
        )
        .arg(
            Arg::new("sequence")
                .short('s')
                .long("sequence")
                .help("Amino acid sequence")
                .value_name("SEQ"),
        )
        .arg(
            Arg::new("samples")
                .short('n')
                .long("samples")
                .help("Number of samples")
                .value_name("NUM")
                .value_parser(clap::value_parser!(u64)),
        )
        .arg(
            Arg::new("forcefield")
                .short('f')
                .long("forcefield")
                .help("Force field (amber03, amber10, amber96, amber99sb)")
                .value_name("FF"),
        )
        .get_matches();

    // Handle config file generation
    if matches.get_one::<bool>("config").is_some_and(|&v| v) {
        generate_config_file()?;
        return Ok(());
    }

    // Try to get parameters in this order: CLI args -> config file -> show help
    let params = if let Some(params) = get_params_from_cli(&matches) {
        params
    } else if let Some(params) = get_params_from_config() {
        params
    } else {
        println!("No parameters provided via CLI or config file.");
        println!("You can:");
        println!("  1. Provide parameters via command line arguments:");
        println!("     dncs -m molecule -s sequence -n samples -f forcefield");
        println!("  2. Generate a config file:");
        println!("     dncs -c");
        println!("  3. Or use an existing dncs.json file");
        return Ok(());
    };

    // Validate forcefield
    let ff = match params.forcefield.as_str() {
        "amber03" => FF::AMBER03,
        "amber10" => FF::AMBER10,
        "amber96" => FF::AMBER96,
        "amber99sb" => FF::AMBER99SB,
        _ => {
            eprintln!(
                "Unsupported forcefield: {}. Must be one of: amber03, amber10, amber96, amber99sb",
                params.forcefield
            );
            std::process::exit(1);
        }
    };

    // Run simulation
    println!("Generating Samples..");
    let mut system = System::new(&params.sequence, ff.init());
    system.init_parameters();
    system.get_dihedral();

    let mut sample = Sampler::new(system);
    sample.sample(params.n_samples);

    // Create result directory
    let result_dir = format!("Result/{}", params.molecule);
    if fs::metadata(&result_dir).is_ok() {
        fs::remove_dir_all(&result_dir)?;
    }
    fs::create_dir_all(&result_dir)?;

    // Write initial samples
    sample.write_sampled_angles(&format!("Result/{}/Sampled.out", params.molecule));
    sample.to_pdb(&format!("Result/{}/Sampled.pdb", params.molecule));

    println!("Executing Minimizer");

    // Initialize and run parallel minimizer
    let mut mini = Minimizer::new(sample);
    mini.minimize().await;
    mini.conformational_sort().await;

    // Write minimized results
    mini.write_sampled_angles(&format!(
        "Result/{}/{}.out",
        params.molecule, params.molecule
    ))
    .await?;
    mini.to_pdb(&format!(
        "Result/{}/{}.pdb",
        params.molecule, params.molecule
    ))
    .await?;

    Ok(())
}

fn get_params_from_cli(matches: &clap::ArgMatches) -> Option<SimulationParams> {
    let molecule = matches.get_one::<String>("molecule")?;
    let sequence = matches.get_one::<String>("sequence")?;
    let n_samples = matches.get_one::<u64>("samples")?.to_owned() as usize;
    let forcefield = matches.get_one::<String>("forcefield")?;

    Some(SimulationParams {
        molecule: molecule.to_string(),
        sequence: sequence.to_string(),
        n_samples,
        forcefield: forcefield.to_string(),
    })
}

fn get_params_from_config() -> Option<SimulationParams> {
    if !Path::new("dncs.json").exists() {
        return None;
    }

    let config_content = match fs::read_to_string("dncs.json") {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Error reading dncs.json: {}", e);
            return None;
        }
    };

    let config: Value = match serde_json::from_str(&config_content) {
        Ok(value) => value,
        Err(e) => {
            eprintln!("Failed to parse JSON: {}", e);
            return None;
        }
    };

    let generate = match config["Generate"].as_object() {
        Some(object) => object,
        None => {
            eprintln!("Failed to get Generate object from config");
            return None;
        }
    };

    Some(SimulationParams {
        molecule: generate["molecule"].as_str().unwrap_or("").to_string(),
        sequence: generate["sequence"].as_str().unwrap_or("").to_string(),
        n_samples: generate["n_samples"].as_u64().unwrap_or(10) as usize,
        forcefield: generate["forcefield"].as_str().unwrap_or("").to_string(),
    })
}

fn generate_config_file() -> Result<(), std::io::Error> {
    let config = json!({
        "Generate": {
            "molecule": "Sample",
            "sequence": "YGGFM",
            "n_samples": 10,
            "forcefield": "amber99sb"
        }
    });

    let config_content = serde_json::to_string_pretty(&config).unwrap();
    fs::write("dncs.json", config_content)?;
    println!("Generated dncs.json file");
    Ok(())
}
