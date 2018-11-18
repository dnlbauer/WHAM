use std::fmt;

// One histogram
#[derive(Debug,Clone)]
pub struct Histogram {
	// total number of data points stored in the histogram
	pub num_points: u32,

	// histogram bins
	pub bins: Vec<f64>
}

impl Histogram {
	pub fn new(num_points: u32, bins: Vec<f64>) -> Histogram {
		Histogram {num_points, bins}
	}
}

// a set of histograms
#[derive(Debug,Clone)]
pub struct Dataset {
	// number of histogram windows (number of simulations)
	pub num_windows: usize,

	// total number of bins
	pub num_bins: usize,

	// number of bins in each dimension
	pub dimens_lengths: Vec<usize>,

	// min values of the histogram in each dimension
	hist_min: Vec<f64>,

	// max values of the histogram in each dimension
	hist_max: Vec<f64>,

	// width of a bin in unit of its dimension
	bin_width: Vec<f64>,

	// value of kT
	pub kT: f64,

	// histogram for each window
	pub histograms: Vec<Histogram>,

	// flag for cyclic reaction coordinates
	pub cyclic: bool,

	// locations of biases
	bias_pos: Vec<f64>,

	// force constants of biases
	bias_fc: Vec<f64>,

	// bias value cache
	bias: Vec<f64>,

	// histogram weight
	pub weights: Vec<f64>,
}

impl Dataset {

	pub fn new(num_bins: usize, dimens_lengths: Vec<usize>, bin_width: Vec<f64>, hist_min: Vec<f64>, hist_max: Vec<f64>, bias_pos: Vec<f64>, bias_fc: Vec<f64>, kT: f64, histograms: Vec<Histogram>, cyclic: bool) -> Dataset {
		let num_windows = histograms.len();
		let bias: Vec<f64> = vec![0.0; num_bins*num_windows];
		let weights = vec![1.0; num_windows];
		let mut ds = Dataset{
			num_windows,
			num_bins,
			dimens_lengths,
			bin_width,
			hist_min,
			hist_max,
			kT,
			histograms,
			cyclic,
			bias_pos,
			bias_fc,
			bias,
			weights
		};
		for window in 0..num_windows {
			for bin in 0..num_bins {
				let ndx = window * num_bins + bin;
				ds.bias[ndx] = ds.calc_bias(bin, window);
			}
		}
		ds

	}

	pub fn new_weighted(ds: Dataset, weights: Vec<f64>) -> Dataset {
		Dataset {
			weights: weights,
			..ds
		}
	}

	pub fn get_weighted_bin_count(&self, bin: usize) -> f64 {
		self.histograms.iter().enumerate().map(|(idx,h)| self.weights[idx]*h.bins[bin]).sum()
	}

	fn expand_index(&self, bin: usize, lengths: &[usize]) -> Vec<usize> {
    	let mut tmp = bin;
    	let mut idx = vec![0; lengths.len()];
    	for dimen in (1..lengths.len()).rev() {
        	let denom = lengths.iter().take(dimen).fold(1, |s,&x| s*x);
			idx[dimen] = tmp / denom;
        	tmp = tmp % denom;
    	}
    	idx[0] = tmp;
    	idx
	}

	// get center x value for a bin
	pub fn get_coords_for_bin(&self, bin: usize) -> Vec<f64> {
		self.expand_index(bin, &self.dimens_lengths).iter().enumerate().map(|(i, dimen_bin)| {
			self.hist_min[i] + self.bin_width[i]*(*dimen_bin as f64 + 0.5)
		}).collect()
	}

	pub fn get_bias(&self, bin: usize, window: usize) -> f64 {
		let ndx = window * self.num_bins + bin;
		self.bias[ndx]
	}

	// Harmonic bias calculation: bias = 0.5*k(dx)^2
	// if cyclic is true, lowest and highest bins are assumed to be
	// neighbors. This returns exp(U/kT) instead of U for better performance.
	fn calc_bias(&self, bin: usize, window: usize) -> f64 {
		let dimens = self.dimens_lengths.len();
		// index of the bias value depends on the window und dimension
		let bias_ndx: Vec<usize> = (0..dimens)
			.map(|dimen| { window * dimens + dimen }).collect();

		// find the N coords, force constants and bias coords
		let coord = self.get_coords_for_bin(bin);
		let bias_fc: Vec<f64> = bias_ndx.iter().map(|ndx| { self.bias_fc[*ndx] }).collect();
		let bias_pos: Vec<f64> = bias_ndx.iter().map(|ndx| { self.bias_pos[*ndx] }).collect();

		let mut bias_sum = 0.0;
		for i in 0..dimens {
			let mut dist = (coord[i] - bias_pos[i]).abs();
			if self.cyclic { // periodic conditions
				let hist_len = self.hist_max[i] - self.hist_min[i];
				if dist > 0.5 * hist_len {
					dist -= hist_len;
				}
			}
			// store exp(U/kT) for better performance
			bias_sum += 0.5 * bias_fc[i] * dist * dist
		}
		let bias_sum = (-bias_sum/self.kT).exp();
		bias_sum
	}
}

impl fmt::Display for Dataset {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		let mut datapoints: u32 = 0;
		for h in &self.histograms {
			datapoints += h.num_points;
		}
		write!(f, "{} windows, {} datapoints", self.num_windows, datapoints)
    }
}

#[cfg(test)]
mod tests {
	use super::*;
	use super::super::k_B;

	macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            assert!(($x-$y).abs() < $d, "{} != {}", $x, $y)
        }
    }

	fn build_hist() -> Histogram {
		Histogram::new(
			22, // num_points
			vec![1.0, 1.0, 3.0, 5.0, 12.0] // bins
		)
	}

	fn build_hist_set() -> Dataset {
		let h = build_hist();
		Dataset::new(
			5, // num bins
			vec![1],
			vec![1.0], // bin width
			vec![0.0], // hist min
			vec![9.0], // hist max
			vec![4.5], // x0
			vec![10.0], // fc
			300.0*k_B, // kT
			vec![h], // hists
			false // cyclic
		)
	}

	#[test]
	fn calc_bias() {
		let ds = build_hist_set(); // k = 10

		// 3th element -> x=3.5, x0=3.5
		assert_delta!(0.134722337796, ds.calc_bias(3, 0), 0.00000001);

		// 8th element -> x=8.5, x0=3.5
		assert_delta!(1.0, ds.calc_bias(4,0), 0.00000001);
		
		// 1st element -> x=0.5, x0=3.5. non-cyclic!
		assert_delta!(0.0, ds.calc_bias(0,0), 0.0000001);
	}

	#[test]
	fn calc_biascyclic() {
		let mut ds = build_hist_set();
		ds.cyclic = true;

		// 7th element -> x=3.5, x0=3.5
		assert_delta!(0.134722337796, ds.calc_bias(3, 0), 0.00000001);

		// 8th element -> x=4.5, x0=3.5
		assert_delta!(1.0, ds.calc_bias(4, 0), 0.00000001);
		

		// 1th element -> x=0.5, x0=3.5
		// cyclic flag makes bin 0 neighboring bin 9, so the distance is actually 2
		assert_delta!(0.0000000000000117769, ds.calc_bias(0, 0), 0.00000001);

		// 2nd element -> x=1.5, x0=3.5
		assert_delta!(0.00000001, ds.calc_bias(1, 0), 0.00000001);
	}

	#[test]
	fn get_x_for_bin() {
		let ds = build_hist_set();
		let expected: Vec<f64> = vec![0,1,2,3,4,5,6,7,8].iter()
				.map(|x| *x as f64 + 0.5).collect(); 
		for i in 0..9 {
			assert_eq!(expected[i], ds.get_coords_for_bin(i)[0]);
		}
	}

	#[test]
	fn get_bin_count() {
		let ds = Dataset::new(
			5, // num bins
			vec![1],
			vec![1.0, 1.0], // bin width
			vec![0.0, 0.0], // hist min
			vec![5.0, 5.0], // hist max
			vec![7.5, 7.5], // x0
			vec![10.0, 10.0], // fc
			300.0*k_B, // kT
			vec![build_hist(), build_hist()], // hists
			false // cyclic
		);
		assert_delta!(2.0, ds.get_weighted_bin_count(0), 0.0000000001);
		assert_delta!(2.0, ds.get_weighted_bin_count(1), 0.0000000001);
		assert_delta!(6.0, ds.get_weighted_bin_count(2), 0.0000000001);
		assert_delta!(10.0, ds.get_weighted_bin_count(3), 0.0000000001);
		assert_delta!(24.0, ds.get_weighted_bin_count(4), 0.0000000001);
	}
}