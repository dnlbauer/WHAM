#![allow(non_snake_case)]

pub mod io;
pub mod histogram;

use std::error::Error;
use std::result::Result;
use histogram::Dataset;
use std::f64;
use std::fmt;
use std::io::prelude::*;

#[allow(non_upper_case_globals)]
static k_B: f64 = 0.0083144621; // kJ/mol*K

// Application config
#[derive(Debug)]
pub struct Config {
	pub metadata_file: String,
	pub hist_min: Vec<f64>,
	pub hist_max: Vec<f64>,
	pub num_bins: Vec<usize>,
	pub dimens: usize,
	pub verbose: bool,
	pub tolerance: f64,
	pub max_iterations: usize,
	pub temperature: f64,
	pub cyclic: bool,
	pub output: String,
}

impl fmt::Display for Config {
	 fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
	 	write!(f, "Metadata={}, hist_min={:?}, hist_max={:?}, bins={:?}\nverbose={}, tolerance={}, iterations={}, temperature={}, cyclic={:?}" , self.metadata_file, self.hist_min,
	 		self.hist_max, self.num_bins, self.verbose, self.tolerance,
	 		self.max_iterations, self.temperature, self.cyclic)
    }
}

// Checks for convergence between two WHAM iterations. WHAM is considered as
// converged if the absolute difference for the calculated bias offset is 
// smaller then a tolerance value for every simulation window.
fn is_converged(old_F: &[f64], new_F: &[f64], tolerance: f64) -> bool {
	!new_F.iter().zip(old_F.iter())
            .map(|x| { (x.0-x.1).abs() })
            .any(|diff| { diff > tolerance })
}

// estimate the probability of a bin of the histogram set based on F values
// This evaluates the first WHAM equation for each bin
fn calc_bin_probability(bin: usize, ds: &Dataset, F: &[f64]) -> f64 {
    let mut denom_sum: f64 = 0.0;
	let mut bin_count: f64 = 0.0;
	for (window, h) in ds.histograms.iter().enumerate() {
		bin_count += h.bins[bin];
		let bias = -ds.kT*ds.calc_bias(bin, window).ln();
		let bias_offset = ((F[window] - bias) / ds.kT).exp();
		denom_sum += (h.num_points as f64) * bias_offset;
	}

	bin_count / denom_sum
}

// estimate the bias offset F of the histogram based on given probabilities
// This evaluates the second WHAM equation for each window
fn calc_window_F(window: usize, ds: &Dataset, P: &[f64]) -> f64 {
    (0..ds.num_bins).zip(P.iter()) // zip bins and P
		.filter_map(|bin_and_prob: (usize, &f64)| {
            if bin_and_prob.1 == &0.0 { // skip zeros for speed
                None
            } else {
                Some(bin_and_prob.1 * ds.calc_bias(bin_and_prob.0, window))
            }
		}).sum()
}

// One full WHAM iteration includes calculation of new probabilities P and
// new bias offsets F based on previous bias offsets F_prev. This updates
// the values in vectors F and P
fn perform_wham_iteration(ds: &Dataset, F_prev: &[f64], F: &mut [f64], P: &mut [f64]) {
	// reset bias offsets
	for window in 0..ds.num_windows {
		F[window] = 0.0;
	}

	// for bin in 0..ds.num_bins {
	// 	let x = get_x_for_bin(bin, ds.hist_min, ds.bin_width);
	// 	let mut num = 0.0;
	// 	let mut denom = 0.0;

	// 	for window in 0..ds.num_windows {
	// 		match ds.histograms[window].get_bin_count(bin) {
	// 			Some(c) => num += c,
	// 			_ => {}
	// 		}
	// 		let bias = calc_bias(
	// 			ds.bias_fc[window],
	// 			ds.bias_pos[window],
	// 			x);
	// 		let bf = ((F_prev[window]-bias) / ds.kT).exp();
	// 		denom += ds.histograms[window].num_points as f64* bf
	// 	}
	// 	P[bin] = num / denom;
	
	// 	for window in 0..ds.num_windows {
	// 		let bias = calc_bias(
	// 				ds.bias_fc[window],
	// 				ds.bias_pos[window],
	// 				x);
	// 		let bf = (-bias/ds.kT).exp() * P[bin];
	// 		F[window] += bf;
	// 	}
	// }

	// for window in 0..ds.num_windows {
	// 	F[window] = -ds.kT * F[window].ln();
	// }

	// evaluate first WHAM equation for each bin to
	// estimage probabilities based on previous offsets (F_prev)
	for bin in 0..ds.num_bins {
		P[bin] = calc_bin_probability(bin, ds, F_prev);
	}

	// evaluate second WHAM equation for each window to
	// estimate new bias offsets from propabilities
	for window in 0..ds.num_windows {
		F[window] = calc_window_F(window, ds, P);
	}

}

// get average difference between two bias offset sets
fn diff_avg(F: &[f64], F_prev: &[f64]) -> f64 {
	let mut F_sum: f64 = 0.0;
	for i in 0..F.len() {
		F_sum += (F[i]-F_prev[i]).abs()
	}
	F_sum / F.len() as f64
} 


// calculate the normalized free energy from normalized probability values
fn free_energy(ds: &Dataset, P: &[f64], A: &mut [f64]) {
	let mut bin_min = f64::MAX;

	// Free energy calculation
	for bin in 0..ds.num_bins {
		A[bin] = -ds.kT*P[bin].ln();
		if A[bin] < bin_min {
			bin_min = A[bin];
		}
	}

	// Make A relative to minimum
	for bin in 0..ds.num_bins {
		A[bin] -= bin_min;
	}
}

pub fn run(cfg: &Config) -> Result<(), Box<Error>>{
	println!("Supplied WHAM options: {}", &cfg);

	// read input data into the histograms object
	println!("Reading input files.");

	let histograms = io::read_data(&cfg) // TODO nicer error handling for this
		.expect("No datapoints in histogram boundaries.");
	println!("{}",&histograms);

	// allocate only once for better performance
	let mut F_prev: Vec<f64> = vec![f64::INFINITY; histograms.num_windows]; 
	let mut F: Vec<f64> = vec![0.0; histograms.num_windows]; 
	let mut P: Vec<f64> = vec![f64::NAN; histograms.num_bins];
	let mut A: Vec<f64> = vec![f64::NAN; histograms.num_bins];

	// perform WHAM until convergence
	let mut iteration = 0;
	while !is_converged(&F_prev, &F, cfg.tolerance) && iteration < cfg.max_iterations {
		iteration += 1;
		// store F values before the next iteration
		F_prev.copy_from_slice(&F[..]);

		// perform wham iteration and update F
		perform_wham_iteration(&histograms, &F_prev, &mut F, &mut P);

		// output some stats during calculation
		if iteration % 10 == 0 {
			println!("Iteration {}: dF={}", &iteration, &diff_avg(&F_prev, &F));
		}

		// Dump free energy and bias offsets
		if iteration % 100 == 0 {
			free_energy(&histograms, &mut P, &mut A);
			dump_state(&histograms, &F, &F_prev, &P, &A);
		}
	}
	
	// Normalize P
	let P_sum: f64 = P.iter().sum();
	P.iter_mut().map(|p| *p /= P_sum).count();

	// final free energy calculation and state dump
	println!("Finished. Dumping final PMF");
	free_energy(&histograms, &mut P, &mut A);
	dump_state(&histograms, &F, &F_prev, &P, &A);

	if iteration == cfg.max_iterations {
		println!("!!!!! WHAM not converged! (max iterations reached) !!!!!");
	}

	io::write_results(&cfg.output, &histograms, &A, &P)?;
	
	Ok(())
}

fn dump_state(ds: &Dataset, F: &[f64], F_prev: &[f64], P: &[f64], A: &[f64]) {
	let out = std::io::stdout();
    let mut lock = out.lock();
	writeln!(lock, "# PMF");
	writeln!(lock, "#x\t\tFree Energy\t\tP(x)");
	for bin in 0..ds.num_bins {
		let x = ds.get_coords_for_bin(bin)[0]; // TODO
		writeln!(lock, "{:9.5}\t{:9.5}\t{:9.5}", x, A[bin], P[bin]);
	}
	writeln!(lock, "# Bias offsets");
	writeln!(lock, "#Window\t\tF\t\tdF");
	for window in 0..ds.num_windows {
		writeln!(lock, "{}\t{:9.5}\t{:8.8}", window, F[window], (F[window]-F_prev[window]).abs());
	}
}


#[cfg(test)]
mod tests {
	use super::histogram::{Dataset,Histogram};
	use std::f64;

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            assert!(($x-$y).abs() < $d, "{} != {}", $x, $y)
        }
    }

	#[test]
	fn is_converged() {
		let new = vec![1.0,1.0];
		let old = vec![0.95, 1.0];
		let tolerance = 0.1;
		let converged = super::is_converged(&old, &new, tolerance);
		assert!(converged);

		let old = vec![0.8, 1.0];
		let converged = super::is_converged(&old, &new, tolerance);
		assert!(!converged);
	}

	fn create_test_ds() -> Dataset {
		let h1 = Histogram::new(10, vec![0.0, 0.0, 3.0, 4.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
		let h2 = Histogram::new(20, vec![0.0, 0.0, 0.0, 3.0, 2.0, 5.0, 10.0, 0.0, 0.0, 0.0, 0.0]);
		Dataset::new(4, vec![1], vec![1.0], vec![0.0], vec![4.0],
                     vec![1.0, 2.0], vec![10.0, 10.0], 2.479, vec![h1, h2], false)
	}

	fn assert_near(a: f64, b: f64, tolerance: f64) {
		let d = (a-b).abs();
		assert!(d <= tolerance, "Values are not close: {}, {}, d={}", &a, &b, &d);
	}

	#[test]
	fn calc_bias_offset() {
		let ds = create_test_ds();
		let probability = vec!(0.959, 0.331, 0.656, 46.750);
        let expected = vec!(0.786289183, 1.10629119);
		for window in 0..ds.num_windows {
			let F = super::calc_window_F(window, &ds, &probability);
            assert_delta!(expected[window], F, 0.0000001);
        }
	}

	#[test]
    #[ignore] // TODO
	fn calc_bin_probability() {
		let ds = create_test_ds();
		let F = vec!(0.0, 0.0);
		let expected = vec!(0.959, 0.331, 0.656, 46.750);
		for b in 0..4 {
			let p = super::calc_bin_probability(b, &ds, &F);
			assert_near(expected[b], p, 0.001);
		}

		let F = vec!(1.0, 1.0);
		let expected = vec!(0.641, 0.221, 0.439, 31.232);
		for b in 0..4 {
			let p = super::calc_bin_probability(b, &ds, &F);
			assert_near(expected[b], p, 0.001);
		}
	}

	#[test]
    #[ignore] // TODO
	fn perform_wham_iteration() {
		let ds = create_test_ds();
		let prev_F = vec![0.0; ds.num_windows];
		let mut F = vec![0.0; ds.num_windows];
		let mut P =  vec![f64::NAN; ds.num_bins];
		super::perform_wham_iteration(&ds, &prev_F, &mut F, &mut P);
		let expected_F = vec!(0.5948, -0.2513);
		let expected_P = vec!(0.959, 0.331, 0.656, 46.750);
		for bin in 0..ds.num_bins {
			assert_near(expected_P[bin], P[bin], 0.01)
		}
		for window in 0..ds.num_windows {
			assert_near(expected_F[window], F[window], 0.01)	
		}
		
	}
}