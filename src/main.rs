/*
 * Digital Nets Conformational Sampling (DNCS)
 * Copyright [2024] [Abraham Rebairo J., Satheeshkumar S, Sam Paul D., Stephen A. and FiveAngstromLLP]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

use clap::{Arg, Command};
use libdncs::*;
use std::sync::Arc;

fn main() {
    let matches = Command::new("DNCS")
        .version("1.0.0")
        .about("Digital Nets Conformational Sample (DNCS)")
        .subcommand(
            Command::new("sample")
            .about("Generate conformational samples")
                .arg(
                    Arg::new("folder")
                        .short('F')
                        .long("folder")
                        .help("Folder name")
                        .value_name("NAME")
                        .required(true),
                )
                .arg(
                    Arg::new("sequence")
                        .short('s')
                        .long("sequence")
                        .help("Amino acid sequence")
                        .value_name("SEQ")
                        .required(true),
                )
                .arg(
                    Arg::new("samples")
                        .short('n')
                        .long("samples")
                        .help("Number of samples")
                        .value_name("NUM")
                        .value_parser(clap::value_parser!(usize))
                        .required(true),
                )
                .arg(
                    Arg::new("forcefield")
                        .short('f')
                        .long("forcefield")
                        .help("Force field (amber03.xml, amber10.xml, amber96.xml, amber99sb.xml, amberfb15.xml)")
                        .value_name("FF")
                        .required(true),
                )
                .arg(
                    Arg::new("grid")
                        .short('g')
                        .long("grid")
                        .help("No of grid to divide sample space")
                        .value_name("grid")
                        .value_parser(clap::value_parser!(usize))
                        .required(true),
                )
        )
        .subcommand(
            Command::new("minimize")
            .about("Minimize generated conformers")
            .arg(
                Arg::new("folder")
                    .short('F')
                    .long("folder")
                    .help("Folder name")
                    .value_name("NAME")
                    .required(true),
            )
            .arg(
                Arg::new("forcefield")
                    .short('f')
                    .long("forcefield")
                    .help("Force field (amber03.xml, amber10.xml, amber96.xml, amber99sb.xml, amberfb15.xml)")
                    .value_name("FF")
                    .required(true),
            )
            .arg(
                Arg::new("minimize")
                    .short('m')
                    .long("minimize")
                    .help("Minimize Top N conformations")
                    .value_name("NUM")
                    .value_parser(clap::value_parser!(usize))
                    .required(true),
            )
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sample", sample_matches)) => {
            let folder = sample_matches
                .get_one::<String>("folder")
                .expect("Folder is required")
                .to_string();
            let sequence = sample_matches
                .get_one::<String>("sequence")
                .expect("Sequence is required")
                .to_string();
            let n_samples = *sample_matches
                .get_one::<usize>("samples")
                .expect("Samples is required");
            let forcefield_str = sample_matches
                .get_one::<String>("forcefield")
                .expect("Forcefield is required")
                .as_str();
            let grid = *sample_matches
                .get_one::<usize>("grid")
                .expect("Grid is required");

            let ff = match FF::from_str(forcefield_str.trim_end_matches(".xml")) {
                Some(ff) => ff,
                None => {
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
                        forcefield_str
                    );
                    std::process::exit(1);
                }
            };

            let mut sys = System::new(&sequence, ff.init());
            sys.init_parameters();
            let mut sample = Sampler::new(Arc::new(sys), grid, folder);
            sample.sample(n_samples);
        }
        Some(("minimize", minimize_matches)) => {
            let folder = minimize_matches
                .get_one::<String>("folder")
                .expect("Folder is required")
                .to_string();
            let forcefield_str = minimize_matches
                .get_one::<String>("forcefield")
                .expect("Forcefield is required")
                .as_str();
            let minimize = *minimize_matches
                .get_one::<usize>("minimize")
                .expect("Minimize is required");

            let ff = match FF::from_str(forcefield_str.trim_end_matches(".xml")) {
                Some(ff) => ff,
                None => {
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
                        forcefield_str
                    );
                    std::process::exit(1);
                }
            };

            let mut minize = Minimizer::new(folder, ff.init(), minimize);
            minize.minimize_all();
        }
        _ => {
            println!("Use -h or --help for more information")
        }
    }
}
