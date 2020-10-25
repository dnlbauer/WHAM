extern crate wham;
#[macro_use]
extern crate clap;
extern crate rand;

use rand::prelude::*;
use clap::App;
use wham::Config;
use wham::errors::*;
use std::process;

// Parse command line arguments into a Config struct
fn cli() -> Result<Config> {
	let yaml = load_yaml!("cli.yml");
	let matches = App::from_yaml(yaml).get_matches();
	let metadata_file = matches.value_of("metadata").unwrap().to_string();
	let verbose: bool = matches.is_present("verbose");
	let temperature: f64 = matches.value_of("temperature").unwrap().parse()
		.chain_err(|| "Cannot read temperature.")?;
	let tolerance: f64 = matches.value_of("tolerance").unwrap_or("0.000001").parse()
		.chain_err(|| "Cannot read tolerance.")?;
	let max_iterations: usize = matches.value_of("iterations").unwrap_or("100000").parse()
		.chain_err(|| "Cannot parse iterations.")?;
	let output = matches.value_of("output").unwrap_or("wham.out").to_string();
    let cyclic: bool = matches.is_present("cyclic");

	let hist_min: Vec<f64> = matches.value_of("min_hist").unwrap()
        .split(',').map(|x| {
            if x.to_ascii_lowercase() == "pi" {
                std::f64::consts::PI
            } else if x.to_ascii_lowercase() == "-pi" {
                -std::f64::consts::PI
            } else {
                x.parse().unwrap()
            }
        }).collect();
	let hist_max: Vec<f64> = matches.value_of("max_hist").unwrap()
        .split(',').map(|x| {
            if x.to_ascii_lowercase() == "pi" {
                std::f64::consts::PI
            } else if x.to_ascii_lowercase() == "-pi" {
                -std::f64::consts::PI
            } else {
                x.parse().unwrap()
            }
        }).collect();
	let num_bins: Vec<usize> = matches.value_of("bins").unwrap()
        .split(',').map(|x| { x.parse().unwrap() }).collect();
	let bootstrap: usize = matches.value_of("bootstrap").unwrap_or("0").parse()
		.chain_err(|| "Cannot parse bootstrap iteration.")?;
    let bootstrap_seed: u64 = matches.value_of("bootstrap_seed")
        .unwrap_or({
            let mut rng = rand::thread_rng();
            &rng.gen::<u32>().to_string()
        }).parse()
        .chain_err(|| "Cannot parse bootstrap iteration.")?;
    let start: f64 = matches.value_of("start").unwrap_or("0").parse()
        .chain_err(|| "Cannot parse start time.")?;
    let end: f64 = matches.value_of("end").unwrap_or("1e+20").parse()
        .chain_err(|| "Cannot parse end time.")?;

    let uncorr: bool = matches.is_present("uncorr");
     
    if num_bins.len() != hist_max.len() || num_bins.len() != hist_max.len() {
        eprintln!("Input dimensions do not match (min: {}, max: {}, bins: {})",
                  hist_min.len(), hist_max.len(), num_bins.len());
        process::exit(1);
    }

    let dimens = num_bins.len();

	Ok(wham::Config{metadata_file, hist_min, hist_max, num_bins, dimens,
		verbose, tolerance, max_iterations, temperature, cyclic, output,
		bootstrap, bootstrap_seed, start, end, uncorr})
}

fn main() {

	let cfg = cli().expect("Failed to parse CLI.");
	if let Err(error) = wham::run(&cfg) {
		eprintln!("Error: {}", error);

		for e in error.iter().skip(1) {
			eprintln!("Reason: {}", e)
		}
		process::exit(1);
	}
}
