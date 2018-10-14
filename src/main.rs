extern crate wham;
#[macro_use]
extern crate clap;

use clap::App;
use wham::Config;
use std::error::Error;
use std::result::Result;
use std::process;

// Parse command line arguments into a Config struct
fn cli() -> Result<Config, Box<Error>> {
	let yaml = load_yaml!("cli.yml");
	let matches = App::from_yaml(yaml).get_matches();
	let metadata_file = matches.value_of("metadata").unwrap().to_string();
	let verbose: bool = matches.is_present("verbose");
	let temperature: f64 = matches.value_of("temperature").unwrap().parse()?;
	let tolerance: f64 = matches.value_of("tolerance").unwrap_or("0.000001").parse()?;
	let max_iterations: usize = matches.value_of("iterations").unwrap_or("100000").parse()?;
	let output = matches.value_of("output").unwrap_or("wham.out").to_string();
    let cyclic: bool = matches.is_present("cyclic");

	let hist_min: Vec<f64> = matches.value_of("min_hist").unwrap()
        .split(',').map(|x| { x.parse().unwrap() }).collect();
	let hist_max: Vec<f64> = matches.value_of("max_hist").unwrap()
        .split(',').map(|x| { x.parse().unwrap() }).collect();
	let num_bins: Vec<usize> = matches.value_of("bins").unwrap()
        .split(',').map(|x| { x.parse().unwrap() }).collect();

    if num_bins.len() != hist_max.len() || num_bins.len() != hist_max.len() {
        eprintln!("Input dimensions do not match (min: {}, max: {}, bins: {})",
                  hist_min.len(), hist_max.len(), num_bins.len());
        process::exit(1);
    }

    let dimens = num_bins.len();

	Ok(wham::Config{metadata_file, hist_min, hist_max, num_bins, dimens,
		verbose, tolerance, max_iterations, temperature, cyclic, output})
}

fn main() {
	let cfg = cli().expect("Failed to parse CLI.");
	match wham::run(&cfg) {
		Err(_) => process::exit(1),
		_ => {}
	}
}
