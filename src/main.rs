// Copyright 2024 Abraham Rebairo J., Satheeshkumar S, Sam Paul D., Stephen A. and FiveAngstromLLP
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use clap::{Arg, Command};
use libdncs::*;
use serde_json::{json, Value};
use std::fs;
use std::path::Path;
use std::sync::Arc;

struct SamplingParams {
    molecule: String,
    sequence: String,
    n_samples: usize,
    forcefield: String,
    minimize: bool,
    grid: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("DNCS")
        .version("1.0.0")
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
                .short('N')
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
                .help("Force field (amber03.xml, amber10.xml, amber96.xml, amber99sb.xml, amberfb15.xml)")
                .value_name("FF"),
        )
        .arg(
            Arg::new("minimize")
                .short('m')
                .long("minimize")
                .help("Bool parameter to ensure minimize")
                .action(clap::ArgAction::SetTrue)
        )
        .arg(
            Arg::new("grid")
                .short('g')
                .long("grid")
                .help("No of grid to divide sample")
                .value_name("grid")
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
        "amber03.xml" => FF::AMBER03,
        "amber10.xml" => FF::AMBER10,
        "amber96.xml" => FF::AMBER96,
        "amber99sb.xml" => FF::AMBER99SB,
        "amberfb15.xml" => FF::AMBERFB15,
        _ => {
            eprintln!(
                "
Unsupported forcefield: {}.
Must be one of below:
    - amber03.xml
    - amber10.xml
    - amber96.xml
    - amber99sb.xml
    - amberfb15.xml
",
                params.forcefield
            );
            std::process::exit(1);
        }
    };

    // Run Sampling
    println!(
        "Generating Best {} Samples from 2048 Samples..",
        &params.n_samples
    );
    let mut system = System::new(&params.sequence, ff.init());
    system.init_parameters();

    let mut sample = Sampler::new(Arc::new(system), params.grid);
    sample.sample(params.n_samples);

    // Create result directory
    let result_dir = format!("Result/{}", params.molecule);
    if fs::metadata(&result_dir).is_ok() {
        fs::remove_dir_all(&result_dir)?;
    }
    fs::create_dir_all(&result_dir)?;

    // Write initial samples
    sample.write_angles(&format!("Result/{}/Sampled.out", params.molecule));
    sample.to_pdb(&format!("Result/{}/Sampled.pdb", params.molecule));

    if params.minimize {
        println!("Executing Minimizer");

        // Initialize and run parallel minimizer
        let mut mini = Minimizer::new(sample);
        mini.minimize();
        mini.conformational_sort();

        // Write minimized results
        mini.write_angles(&format!(
            "Result/{}/{}.out",
            params.molecule, params.molecule
        ))?;
        mini.to_pdb(&format!(
            "Result/{}/{}.pdb",
            params.molecule, params.molecule
        ))?;
    }
    Ok(())
}

fn get_params_from_cli(matches: &clap::ArgMatches) -> Option<SamplingParams> {
    let molecule = matches.get_one::<String>("molecule")?;
    let sequence = matches.get_one::<String>("sequence")?;
    let n_samples = matches.get_one::<u64>("samples")?.to_owned() as usize;
    let grid = matches.get_one::<u64>("grid")?.to_owned() as usize;
    let forcefield = matches.get_one::<String>("forcefield")?;
    let minimize = matches.get_one::<bool>("minimize");

    Some(SamplingParams {
        molecule: molecule.to_string(),
        sequence: sequence.to_string(),
        n_samples,
        grid,
        forcefield: forcefield.to_string(),
        minimize: minimize.unwrap_or(&false).to_owned(),
    })
}

fn get_params_from_config() -> Option<SamplingParams> {
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

    Some(SamplingParams {
        molecule: generate["molecule"].as_str().unwrap_or("").to_string(),
        sequence: generate["sequence"].as_str().unwrap_or("").to_string(),
        n_samples: generate["n_samples"].as_u64().unwrap_or(10) as usize,
        grid: generate["grid"].as_u64().unwrap_or(10) as usize,
        forcefield: generate["forcefield"].as_str().unwrap_or("").to_string(),
        minimize: generate["minimize"].as_bool().unwrap_or(false) as bool,
    })
}

fn generate_config_file() -> Result<(), std::io::Error> {
    let config = json!({
        "Generate": {
            "molecule": "Sample",
            "sequence": "YGGFM",
            "n_samples": 10,
            "forcefield": "amberfb15.xml",
            "minimize": true,
            "grid": 4
        }
    });

    let config_content = serde_json::to_string_pretty(&config).unwrap();
    fs::write("dncs.json", config_content)?;
    println!("Generated dncs.json file");
    Ok(())
}
