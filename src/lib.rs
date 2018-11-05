#![allow(non_snake_case)]

#[macro_use]
extern crate error_chain;

pub mod io;
pub mod histogram;

use histogram::Dataset;
use std::f64;
use std::fmt;
use std::io::prelude::*;

// init error chain
pub mod errors { error_chain!{} }
use errors::*;

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
         write!(f, "Metadata={}, hist_min={:?}, hist_max={:?}, bins={:?} verbose={}, tolerance={}, iterations={}, temperature={}, cyclic={:?}", self.metadata_file, self.hist_min, self.hist_max, self.num_bins,
                self.verbose, self.tolerance, self.max_iterations, self.temperature, self.cyclic)
    }
}

// Checks for convergence between two WHAM iterations. WHAM is considered as
// converged if the maximal difference for the calculated bias offsets is
// smaller then a tolerance value.
fn is_converged(old_F: &[f64], new_F: &[f64], tolerance: f64) -> bool {
	!new_F.iter().zip(old_F.iter())
            .map(|x| { (x.0-x.1).abs() })
            .any(|diff| { diff > tolerance })
}

// estimate the probability of a bin of the histogram set based on given bias offsets (F)
// This evaluates the first WHAM equation for each bin.
fn calc_bin_probability(bin: usize, dataset: &Dataset, F: &[f64]) -> f64 {
    let mut denom_sum: f64 = 0.0;
    for (window, h) in dataset.histograms.iter().enumerate() {
		let bias = dataset.calc_bias(bin, window);
        denom_sum += (h.num_points as f64) * bias * F[window];
	}
    dataset.bin_count[bin] / denom_sum
}

// estimate the bias offset F of the histogram based on given probabilities.
// This evaluates the second WHAM equation for each window and returns exp(F/kT)
fn calc_window_F(window: usize, dataset: &Dataset, P: &[f64]) -> f64 {
    let f: f64 = (0..dataset.num_bins).zip(P.iter()) // zip bins and P
        .map(|bin_and_prob: (usize, &f64)| {
            let bias = dataset.calc_bias(bin_and_prob.0, window);
            bin_and_prob.1 * bias
        }).sum();
    1.0/f
}

// One full WHAM iteration includes calculation of new probabilities P and
// new bias offsets F based on previous bias offsets F_prev. This updates
// the values in vectors F and P
fn perform_wham_iteration(dataset: &Dataset, F_prev: &[f64], F: &mut [f64], P: &mut [f64]) {
	// evaluate first WHAM equation for each bin to
	// estimage probabilities based on previous offsets (F_prev)
	for bin in 0..dataset.num_bins {
		P[bin] = calc_bin_probability(bin, dataset, F_prev);
	}

	// evaluate second WHAM equation for each window to
	// estimate new bias offsets from propabilities
	for window in 0..dataset.num_windows {
		F[window] = calc_window_F(window, dataset, P);
	}
}

pub fn perform_wham(cfg: &Config, dataset: &Dataset) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
	// allocate required vectors.
    let mut P: Vec<f64> = vec![f64::NAN; dataset.num_bins]; // bin probability
    let mut F: Vec<f64> = vec![1.0; dataset.num_windows]; // bias offset exp(F/kT)
    let mut F_prev: Vec<f64> = vec![f64::NAN; dataset.num_windows]; // previous bias offset
    let mut F_tmp: Vec<f64> = vec![f64::NAN; dataset.num_windows]; // temp storage for F

    let mut iteration = 0;
    let mut converged = false;

    // perform WHAM until convergence
    while !converged && iteration < cfg.max_iterations {
        iteration += 1;

        // store F values before the next iteration
        F_prev.copy_from_slice(&F);

        // perform wham iteration (this updates F and P)
        perform_wham_iteration(&dataset, &F_prev, &mut F, &mut P);

        // convergence check
        if iteration % 10 == 0 {
            // This backups exp(F/kT) in a temporary vector and calculates true F and F_prev for
            // convergence. Finally, F is restored. F_prev does not need to be restored because
            // its overwritten for the next iteration.
            F_tmp.copy_from_slice(&F);
            for f in F.iter_mut() { *f = -dataset.kT * f.ln() }
            for f in F_prev.iter_mut() { *f = -dataset.kT * f.ln() }
            converged = is_converged(&F_prev, &F, cfg.tolerance);

            println!("Iteration {}: dF={}", &iteration, &diff_avg(&F_prev, &F));
            F.copy_from_slice(&F_tmp);
        }

        // Dump free energy and bias offsets
        //if iteration % 100 == 0 {
        //   free_energy(&histograms, &mut P, &mut A);
        //    dump_state(&histograms, &F, &F_prev, &P, &A);
        //}
    }

    // Normalize P to sum(P) = 1.0
    let P_sum: f64 = P.iter().sum();
    P.iter_mut().map(|p| *p /= P_sum).count();

    if iteration == cfg.max_iterations {
		bail!("WHAM not converged! (max iterations reached)");
    }

	Ok((P, F, F_prev))
}

pub fn run(cfg: &Config) -> Result<()>{
    println!("Supplied WHAM options: {}", &cfg);

    println!("Reading input files.");
    // TODO Better error handling with nice error messages instead of a panic!
    let dataset = io::read_data(&cfg).chain_err(|| "Failed to create histogram.")?;
    println!("{}", &dataset);

    let (P, F, F_prev) = perform_wham(&cfg, &dataset)?;

    // calculate free energy and dump state
    println!("Finished. Dumping final PMF");
	let free_energy = calc_free_energy(&dataset, &P);
    dump_state(&dataset, &F, &F_prev, &P, &free_energy);

    io::write_results(&cfg.output, &dataset, &free_energy, &P)
		.chain_err(|| "Could not write results to output file")?;

    Ok(())
}


// get average difference between two bias offset sets
fn diff_avg(F: &[f64], F_prev: &[f64]) -> f64 {
	let mut F_sum: f64 = 0.0;
	for i in 0..F.len() {
		F_sum += (F[i]-F_prev[i]).abs()
	}
	F_sum / F.len() as f64
}

// calculate the normalized free energy from probability values
fn calc_free_energy(dataset: &Dataset, P: &[f64]) -> Vec<f64> {
    let mut minimum = f64::MAX;
	let mut free_energy: Vec<f64> = P.iter()
        .map(|p| {
            -dataset.kT * p.ln()
        })
        .inspect(|free_e| {
            if free_e < &minimum {
                minimum = *free_e;
            }
        })
        .collect();

    for e in free_energy.iter_mut() {
        *e -= minimum
    }
    free_energy
}

fn dump_state(dataset: &Dataset, F: &[f64], F_prev: &[f64], P: &[f64], A: &[f64]) {
	let out = std::io::stdout();
    let mut lock = out.lock();
	writeln!(lock, "# PMF");
	writeln!(lock, "#bin\t\tFree Energy\t\tP(x)");
	for bin in 0..dataset.num_bins {
		writeln!(lock, "{:9.5}\t{:9.5}\t{:9.5}", bin, A[bin], P[bin]);
	}
	writeln!(lock, "# Bias offsets");
	writeln!(lock, "#Window\t\tF\t\tdF");
	for window in 0..dataset.num_windows {
		writeln!(lock, "{}\t{:9.5}\t{:8.8}", window, F[window], (F[window]-F_prev[window]).abs());
	}
}


#[cfg(test)]
mod tests {
	use super::histogram::{Dataset,Histogram};
	use std::f64;
    use super::k_B;

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            assert!(($x-$y).abs() < $d, "{} != {}", $x, $y)
        }
    }


	fn create_test_dataset() -> Dataset {
		let h1 = Histogram::new(10, vec![0.0, 1.0, 1.0, 8.0, 0.0]);
		let h2 = Histogram::new(10, vec![0.0, 0.0, 8.0, 1.0, 1.0]);
		Dataset::new(5, vec![5], vec![1.0], vec![0.0], vec![4.0],
                     vec![1.0, 1.0], vec![10.0, 10.0], 300.0*k_B, vec![h1, h2], false)
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

    #[test]
	fn calc_bin_probability() {
		let dataset = create_test_dataset();
		let F = vec![1.0; dataset.num_bins]  ;
        let expected = vec!(0.0, 0.0825296687031316, 40.92355847097493,
                            124226.70003377, 2308526035.5283747);
		for b in 0..dataset.num_bins {
			let p = super::calc_bin_probability(b, &dataset, &F);
			assert_delta!(expected[b], p, 0.0000001);
		}
	}

    #[test]
	fn calc_bias_offset() {
		let dataset = create_test_dataset();
		let probability = vec!(0.0, 0.1, 0.2, 0.3, 0.4);
        let expected = vec!(15.927477169990633, 15.927477169990633);
		for window in 0..dataset.num_windows {
			let F = super::calc_window_F(window, &dataset, &probability);
            assert_delta!(expected[window], F, 0.0000001);
        }
	}

	#[test]
	fn perform_wham_iteration() {
		let dataset = create_test_dataset();
		let prev_F = vec![1.0; dataset.num_windows];
		let mut F = vec![f64::NAN; dataset.num_windows];
		let mut P =  vec![f64::NAN; dataset.num_bins];
		super::perform_wham_iteration(&dataset, &prev_F, &mut F, &mut P);
        let expected_F = vec!(1.0, 1.0);
		let expected_P = vec!(0.0, 0.0825296687031316, 40.92355847097493,
                            124226.70003377, 2308526035.5283747);
		for bin in 0..dataset.num_bins {
			assert_delta!(expected_P[bin], P[bin], 0.01)
		}
		for window in 0..dataset.num_windows {
			assert_delta!(expected_F[window], F[window], 0.01)	
		}
		
	}
}