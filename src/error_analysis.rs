use rand;
use super::histogram::{Dataset};
use super::perform_wham;
use super::Config;
use rgsl::statistics;

// returns a set of num_windows continious weights by
// a) generate num_windows-1 random variables and sort them
// b) each weight n is the difference between n+1 and n, where n0=0 and nN+1=1
fn generate_random_weights(num_windows: usize) -> Vec<f64> {
    // create a list of num_windows - 1 sorted random numbers and append/prepend 0 and 1
    let mut tmp = (0..num_windows-1).map(|_| rand::random::<f64>()).collect::<Vec<f64>>();
    tmp.sort_by(|a,b| { a.partial_cmp(b).unwrap() });
    let mut rnds = vec![0.0];
    rnds.append(&mut tmp);
    rnds.append(&mut vec![1.0]);

    // weights of window i is the difference between rnd[i+1] and rnd[i]
    let mut weights = vec![0.0; num_windows];
    for i in 0..num_windows {
        weights[i] = rnds[i+1] - rnds[i]
    }
    return weights
}

// Generate a random weighted dataset from the given dataset by changing the weights
fn generate_random_weighted_dataset(ds: Dataset) -> Dataset {
    let weights = generate_random_weights(ds.num_windows);
    Dataset::new_weighted(ds, weights)
}

// Perform bootstrap error analysis. This runs the WHAM analysis num_runs times on random weighted
// datasets. The standard deviation is calculated on the bootstrapped probabilities of each bin. The
// standard deviation of the free eneergy is then deduced by error propagation (A_std = kT*1/P*P_std)
pub fn run_bootstrap(cfg: &Config, ds: Dataset, P: &[f64], num_runs: usize) -> (Vec<f64>,Vec<f64>) {
    // Calculate bootstrapped probabilities
    let bootstrapped_Ps: Vec<Vec<f64>> = (0..num_runs).map(|x| {
        println!("Bootstrap run {}/{}", x, num_runs);
        let rnd_weighted_dataset = generate_random_weighted_dataset(ds.clone());
        perform_wham(cfg, &rnd_weighted_dataset).unwrap().0
    }).collect();

    // Evaulate standard deviation of P per bin
    let mut P_std = vec![0.0; ds.num_bins];
    for bin in 0..ds.num_bins {
        let Ps = bootstrapped_Ps.iter().map(|window| window[bin]).collect::<Vec<f64>>();
        P_std[bin] = statistics::sd(&Ps, 1, num_runs);
    }

    // A_std by error propagation
    let A_std = P_std.iter().zip(P.iter()).map(|(std,P)| ds.kT*1.0/P*std).collect();
    (P_std, A_std)
}

#[cfg(tests)]
mod tests {
    use super::*;
    use super::super::k_B;
    use super::super::histogram::Histogram;

    fn build_hist() -> Histogram {
		Histogram::new(
			22, // num_points
			vec![1.0, 1.0, 3.0, 5.0, 12.0] // bins
		)
	}

	fn build_hist_set() -> Dataset {
		let h1 = build_hist();
		let h2 = build_hist();
		let h3 = build_hist();
		Dataset::new(
			5, // num bins
			vec![3],
			vec![1.0], // bin width
			vec![0.0], // hist min
			vec![9.0], // hist max
			vec![4.5, 4.5, 4.5], // x0
			vec![10.0, 10.0, 10.0], // fc
			300.0*k_B, // kT
			vec![h1, h2, h3], // hists
			false // cyclic
		)
	}

    #[test]
    fn random_weights() {
        let num_windows = 5;
        let weights = generate_random_weights(num_windows);
        assert_eq!(num_windows, weights.len());
        for w in weights {
            assert!(0.0 < w && w < 1.0);
        }
    }

    #[test]
    fn random_weighted_dataset() {
        let ds = build_hist_set();
        let rnd_weights_ds = generate_random_weighted_dataset(ds);
        println!("{:?}", rnd_weights_ds.weights);
        for w in rnd_weights_ds.weights {
            assert!(w > 0.0);
            assert!(w < 1.0);
        }
    }

}
