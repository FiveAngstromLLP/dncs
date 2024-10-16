mod forcefield;
mod parser;
mod sampling;
mod system;

use crate::sampling::Sampler;
use crate::system::System;
use clap::{Arg, Command};
use std::fs;
use std::io::Write;
use toml;
use zip;
use zip::write::SimpleFileOptions;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("DNCS")
        .version("0.1.3")
        .about("Digital Nets Conformational Sample (DNCS)")
        .arg(
            Arg::new("config")
                .short('c')
                .long("config")
                .help("Generates dncs.toml file")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    if matches.get_one::<bool>("config").is_some_and(|&v| v) {
        generate_config_file()?;
        return Ok(());
    }

    // Read config from TOML file
    let config_content = match fs::read_to_string("dncs.toml") {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Error reading dncs.toml: {}", e);
            std::process::exit(1);
        }
    };
    let config: toml::Value = match toml::from_str(&config_content) {
        Ok(value) => value,
        Err(e) => {
            eprintln!("Failed to parse TOML: {}", e);
            std::process::exit(1);
        }
    };

    let generate = match config["Generate"].as_table() {
        Some(table) => table,
        None => {
            eprintln!("Failed to get Generate table from config");
            std::process::exit(1);
        }
    };
    let molecule = generate["molecule"].as_str().unwrap_or("Sample");
    let sequence = generate["sequence"].as_str().unwrap_or("YGGFM");
    let include_sidechain = generate["include_sidechain"].as_bool().unwrap_or(false);
    let samples = generate["n_samples"].as_integer().unwrap_or(10) as usize;

    let mut system = System::new(sequence);
    system.init_parameters();
    system.get_dihedralatoms(include_sidechain);
    let mut sample = Sampler::new(system);
    sample.sample(samples);
    sample.conformational_sort();

    match fs::create_dir_all(&format!("Result/{}", molecule)) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to create result directory: {}", e);
            std::process::exit(1);
        }
    }

    sample.write_sampled_angles(&format!("Result/{}/angles.out", molecule));
    sample.to_pdb(&format!("Result/{}/{}.pdb", molecule, molecule));

    // Create a zip file
    let file = match fs::File::create(&format!("Result/{}/result.zip", molecule)) {
        Ok(file) => file,
        Err(e) => {
            eprintln!("Failed to create zip file: {}", e);
            std::process::exit(1);
        }
    };
    let mut zip = zip::ZipWriter::new(file);

    // Define file options
    let options = SimpleFileOptions::default().compression_method(zip::CompressionMethod::Stored);

    // Add angles.out to the zip file
    let angles_content = match fs::read(&format!("Result/{}/angles.out", molecule)) {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Failed to read angles.out: {}", e);
            std::process::exit(1);
        }
    };
    match zip.start_file("angles.out", options) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to start angles.out file in zip: {}", e);
            std::process::exit(1);
        }
    }
    match zip.write_all(&angles_content) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to write angles.out to zip: {}", e);
            std::process::exit(1);
        }
    }

    // Add molecule.pdb to the zip file
    let pdb_content = match fs::read(&format!("Result/{}/{}.pdb", molecule, molecule)) {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Failed to read {}.pdb: {}", molecule, e);
            std::process::exit(1);
        }
    };
    match zip.start_file(format!("{}.pdb", molecule), options) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to start {}.pdb file in zip: {}", molecule, e);
            std::process::exit(1);
        }
    }
    match zip.write_all(&pdb_content) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to write {}.pdb to zip: {}", molecule, e);
            std::process::exit(1);
        }
    }

    // Finish creating the zip file
    match zip.finish() {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to finish zip file: {}", e);
            std::process::exit(1);
        }
    }
    println!("Done");
    Ok(())
}

fn generate_config_file() -> Result<(), std::io::Error> {
    let config_content = r#"
[Generate]
molecule = "Sample"
sequence = "AAAAA"
include_sidechain = true
n_samples = 10
"#;
    fs::write("dncs.toml", config_content)?;
    println!("Generated dncs.toml file");
    Ok(())
}
