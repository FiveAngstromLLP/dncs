use clap::{Parser, Subcommand};
use libdncs::*;
use std::sync::Arc;

mod tui;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    subcommand: SubCommands,
}

#[derive(Subcommand, Debug)]
enum SubCommands {
    #[command(about = "Generate conformational samples")]
    Sample(SampleArgs),
    #[command(about = "Minimize generated conformers")]
    Minimize(MinimizeArgs),
    #[command(about = "Generate dihedral angles")]
    Pdbtoangle {
        /// PDB file path
        file: String,
    },
    #[command(about = "TUI for Configurations")]
    Tui,
}

#[derive(Parser, Debug)]
struct SampleArgs {
    /// Folder name
    #[arg(short = 'F', long = "folder", value_name = "NAME")]
    folder: String,

    /// Amino acid sequence
    #[arg(short = 's', long = "sequence", value_name = "SEQ")]
    sequence: String,

    /// Number of samples
    #[arg(short = 'n', long = "samples", value_name = "NUM")]
    samples: usize,

    /// Force Field
    #[arg(short = 'f', long = "forcefield", value_name = "FF")]
    forcefield: FF,

    /// Methods
    #[arg(short = 'm', long = "method", value_name = "METHOD")]
    method: Method,

    /// Temperature
    #[arg(short = 't', long = "temperature", value_name = "TEMP")]
    temperature: f64,
}

#[derive(Parser, Debug)]
struct MinimizeArgs {
    /// Folder name
    #[arg(short = 'F', long = "folder", value_name = "NAME")]
    folder: String,

    /// Force Field
    #[arg(short = 'f', long = "forcefield", value_name = "FF")]
    forcefield: FF,

    /// Minimize Top N conformations
    #[arg(short = 'm', long = "minimize", value_name = "NUM")]
    minimize: usize,
}

fn main() {
    let args = Args::parse();

    match args.subcommand {
        SubCommands::Sample(sample_args) => {
            let mut sys = System::new(&sample_args.sequence, sample_args.forcefield.init());
            sys.init_parameters();
            let mut sample = Sampler::new(Arc::new(sys), sample_args.method, sample_args.folder);
            sample.sample(sample_args.samples);
        }
        SubCommands::Minimize(minimize_args) => {
            let mut minize = Minimizer::new(
                minimize_args.folder,
                minimize_args.forcefield.init(),
                minimize_args.minimize,
            );
            minize.minimize_all();
        }
        SubCommands::Pdbtoangle { file } => {
            let dihedral_angle = RotateAtDihedral::from_pdb(&file)
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(",");
            println!("{}", dihedral_angle);
        }
        SubCommands::Tui => {
            if let Err(e) = tui::run_peptide_selector() {
                eprintln!("TUI Error: {}", e);
                std::process::exit(1);
            }
        }
    }
}
