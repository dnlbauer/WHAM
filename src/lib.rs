#![allow(non_snake_case)]

#[macro_use]
extern crate error_chain;
extern crate rand;
extern crate rgsl;
extern crate rayon;

pub mod io;
pub mod histogram;
pub mod error_analysis;

use histogram::Dataset;
use std::f64;
use std::fmt;
use std::io::prelude::*;
use rayon::prelude::*;

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
	pub bootstrap: usize,
    pub bootstrap_seed: u64,
    pub start: f64,
    pub end: f64,
}

impl fmt::Display for Config {
	 fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
         write!(f, "Metadata={}, hist_min={:?}, hist_max={:?}, bins={:?}, 
            verbose={}, tolerance={}, iterations={}, temperature={},
            cyclic={:?}, bootstrap={:?}, seed={:?}",
            self.metadata_file, self.hist_min, self.hist_max, self.num_bins,
            self.verbose, self.tolerance, self.max_iterations, self.temperature,
            self.cyclic, self.bootstrap, self.bootstrap_seed)
    }
}

// Checks for convergence between two WHAM iterations. WHAM is considered as
// converged if the maximal difference for the calculated bias offsets is
// smaller then a tolerance value.
fn is_converged(old_F: &[f64], new_F: &[f64], tolerance: f64) -> bool {
    // calculates abs diff between every old and new F and checks if any
    // is larger than tolerance 
	!new_F.iter()
        .zip(old_F.iter())
        .map(|x| { (x.0-x.1).abs() })
        .any(|diff| { diff > tolerance })
}

// estimate the probability of a bin of the histogram set based on given bias
// offsets (F). This evaluates the first WHAM equation for each bin:
// P(x) = \frac {\sum_{i=1}^N{n_i(x)}}
//              {\sum_{i=1}^N{  N_i exp(\beta [F_i - U_{bias,i}(x)])}}
fn calc_bin_probability(bin: usize, dataset: &Dataset, F: &[f64]) -> f64 {
	let mut denom_sum: f64 = 0.0;
	let bin_count: f64 = dataset.get_weighted_bin_count(bin);
    for (window, h) in dataset.histograms.iter().enumerate() {
		let bias = dataset.get_bias(bin, window);
        denom_sum += (dataset.weights[window] * h.num_points as f64)
                    * bias * F[window];
	}
    bin_count / denom_sum
}

// estimate the bias offset F of the histogram based on given probabilities.
// This evaluates the second WHAM equation for each window and returns exp(F/kT).
// exp(F/kT) is not required in intermediate steps so we save some time by not
// calculating it for every iteration. 
// F_i = - 1/\beta ln[\sum_{X_{bins}}{P(x)exp(-\beta U_{bias,i}(x))}]
fn calc_window_F(window: usize, dataset: &Dataset, P: &[f64]) -> f64 {
    let f: f64 = (0..dataset.num_bins).zip(P.iter()) // zip bins and P
        .map(|bin_and_prob: (usize, &f64)| {
            let bias = dataset.get_bias(bin_and_prob.0, window);
            bin_and_prob.1 * bias
        }).sum();
    1.0/f
}

// One full WHAM iteration: calculation of new probabilities P and new bias
// offsets F based on previous bias offsets F_prev. This updates the values in
// vectors F and P.
fn perform_wham_iteration(dataset: &Dataset, F_prev: &[f64], F: &mut Vec<f64>, P: &mut Vec<f64>) {
	// Update P
    // evaluate first WHAM equation for each bin to
	// estimate probabilities based on previous offsets (F_prev))
    (0..dataset.num_bins).into_par_iter()
		.map(|bin| { calc_bin_probability(bin, dataset, F_prev) })
		.collect_into_vec(P);

    // Update F
	// evaluate second WHAM equation for each window to
	// estimate new bias offsets from propabilities
	(0..dataset.num_windows).into_par_iter()
		.map(|window| {calc_window_F(window, dataset, P)} )
		.collect_into_vec(F);
}

// Full WHAM calculation. Calls `perform_wham_iteration` until convergence
// criteria are met or max iterations reached.
pub fn perform_wham(cfg: &Config, dataset: &Dataset)
        -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
	// allocate required vectors.

    // bin probability
    let mut P: Vec<f64> = vec![f64::NAN; dataset.num_bins];
    // bias offset exp(F/kT)
    let mut F: Vec<f64> = vec![1.0; dataset.num_windows];
    // previous bias offset
    let mut F_prev: Vec<f64> = vec![f64::NAN; dataset.num_windows];
    // temp storage for F
    let mut F_tmp: Vec<f64> = vec![f64::NAN; dataset.num_windows];

    let mut iteration = 0;
    let mut converged = false;

    // perform WHAM until convergence
    while !converged && iteration < cfg.max_iterations {
        iteration += 1;

        // store F values before the next iteration
        F_prev.copy_from_slice(&F);

        // perform wham iteration (this updates F and P).
        perform_wham_iteration(&dataset, &F_prev, &mut F, &mut P);

        // convergence check
        if iteration % 10 == 0 {
            // This backups exp(F/kT) in a temporary vector and calculates
            // true F and F_prev for convergence. Finally, F is restored.
            // F_prev does not need to be restored because its overwritten
            // for the next iteration.
            F_tmp.copy_from_slice(&F);
            for f in F.iter_mut() { *f = -dataset.kT * f.ln() }
            for f in F_prev.iter_mut() { *f = -dataset.kT * f.ln() }
            converged = is_converged(&F_prev, &F, cfg.tolerance);

            println!("Iteration {}: dF={}", &iteration, &diff_avg(&F_prev, &F));
            F.copy_from_slice(&F_tmp);
        }
    }

    // Normalize P to sum(P) = 1.0
    let P_sum: f64 = P.iter().sum();
    for p in P.iter_mut() {
        *p /= P_sum; 
    }

    if iteration == cfg.max_iterations {
		bail!("WHAM not converged! (max iterations reached)");
    }

	Ok((P, F, F_prev))
}

pub fn run(cfg: &Config) -> Result<()>{
    println!("Supplied WHAM options: {}", &cfg);

    println!("Reading input files.");
    let dataset = io::read_data(&cfg).chain_err(|| "Failed to create histogram.")?;
    println!("{}", &dataset);

    let (P, F, F_prev) = perform_wham(&cfg, &dataset)?;

	let P_std: Vec<f64>;
	let free_energy_std: Vec<f64>;
	if cfg.bootstrap > 0 {
		let error_est = error_analysis::run_bootstrap(&cfg, dataset.clone(), &P, cfg.bootstrap);
		P_std = error_est.0;
		free_energy_std = error_est.1;
	} else {
		P_std = vec![0.0; P.len()];
		free_energy_std = vec![0.0; P.len()];
	}

    // calculate free energy and dump state
    println!("Finished. Dumping final PMF");
	let free_energy = calc_free_energy(&dataset, &P);
    dump_state(&dataset, &F, &F_prev, &P, &P_std, &free_energy, &free_energy_std);

    io::write_results(&cfg.output, &dataset, &free_energy, &free_energy_std, &P, &P_std)
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

// Print the current WHAM iteration state. Dumps the PMF and associated vectors 
fn dump_state(dataset: &Dataset, F: &[f64], F_prev: &[f64], P: &[f64],
    P_std: &[f64], A: &[f64], A_std: &[f64]) {
	// TODO fix output of F/F_prev
	let out = std::io::stdout();
    let mut lock = out.lock();
	writeln!(lock, "# PMF").unwrap();
	writeln!(lock, "#bin\t\tFree Energy\t\t+/-\t\tP(x)\t\t+/-").unwrap();
	for bin in 0..dataset.num_bins {
		writeln!(lock, "{:9.5}\t{:9.5}\t{:9.5}\t{:9.5}\t{:9.5}",
            bin, A[bin], A_std[bin], P[bin], P_std[bin]).unwrap();
	}
	writeln!(lock, "# Bias offsets").unwrap();
	writeln!(lock, "#Window\t\tF\t\tF_prev").unwrap();
	for window in 0..dataset.num_windows {
		writeln!(lock, "{}\t{:9.5}\t{:8.8}",
            window, F[window], (F[window]-F_prev[window]).abs()).unwrap();
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
        let expected = vec!(0.0, 0.082_529_668_703_131_6, 40.923_558_470_974_93,
                            124_226.700_033_77, 2_308_526_035.528_374_7);
		for b in 0..dataset.num_bins {
			let p = super::calc_bin_probability(b, &dataset, &F);
			assert_delta!(expected[b], p, 0.000_000_1);
		}
	}

    #[test]
	fn calc_bias_offset() {
		let dataset = create_test_dataset();
		let probability = vec!(0.0, 0.1, 0.2, 0.3, 0.4);
        let expected = vec!(15.927_477_169_990_633, 15.927_477_169_990_633);
		for window in 0..dataset.num_windows {
			let F = super::calc_window_F(window, &dataset, &probability);
            assert_delta!(expected[window], F, 0.000_000_1);
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
		let expected_P = vec!(0.0, 0.082_529_668_703_131_6, 40.923_558_470_974_93,
                            124_226.700_033_77, 2_308_526_035.528_374_7);
		for bin in 0..dataset.num_bins {
			assert_delta!(expected_P[bin], P[bin], 0.01)
		}
		for window in 0..dataset.num_windows {
			assert_delta!(expected_F[window], F[window], 0.01)	
		}
		
	}
}