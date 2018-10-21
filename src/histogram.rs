use std::fmt;
use std::cell::RefCell;

// One histogram
#[derive(Debug)]
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
#[derive(Debug)]
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
	bias: RefCell<Vec<Option<f64>>>
}

impl Dataset {

	pub fn new(num_bins: usize, dimens_lengths: Vec<usize>, bin_width: Vec<f64>, hist_min: Vec<f64>, hist_max: Vec<f64>, bias_pos: Vec<f64>, bias_fc: Vec<f64>, kT: f64, histograms: Vec<Histogram>, cyclic: bool) -> Dataset {
		let num_windows = histograms.len();
		let bias: RefCell<Vec<Option<f64>>> = RefCell::new(vec![None; num_bins*num_windows]);
		Dataset{
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
		}
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

	// Harmonic bias calculation: bias = 0.5*k(dx)^2
	// if cyclic is true, lowest and highest bins are assumed to be
	// neighbors
	pub fn calc_bias(&self, bin: usize, window: usize) -> f64 {
		let ndx = window * self.num_bins + bin;
		let mut cache = self.bias.borrow_mut();
		match cache[ndx] {
			Some(val) => val,
			None => {
				// TODO optimize this part!
				let dimens = self.hist_min.len();

				let bias_ndx: Vec<usize> = (0..dimens)
					.map(|dimen| { window * dimens + dimen }).collect();
				let coord = self.get_coords_for_bin(bin);
				let bias_fc: Vec<f64> = bias_ndx.iter().map(|ndx| { self.bias_fc[*ndx] }).collect();
				let bias_pos: Vec<f64> = bias_ndx.iter().map(|ndx| { self.bias_pos[*ndx] }).collect();

				let mut bias_sum = 0.0;
				for i in 0..dimens {
					let mut dist = (coord[i] - bias_pos[i]).abs();
					if self.cyclic {
						let hist_len = self.hist_max[i] - self.hist_min[i];
						if dist > 0.5 * hist_len {
							dist -= hist_len;
						}
					}
					bias_sum += 0.5 * bias_fc[i] * dist * dist
				}
				let bias_sum = (-bias_sum/self.kT).exp();
				cache[ndx] = Some(bias_sum);
				bias_sum
			}
		}
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
			9, // num bins
			vec![1],
			vec![1.0], // bin width
			vec![0.0], // hist min
			vec![9.0], // hist max
			vec![7.5], // x0
			vec![10.0], // fc
			300.0*k_B, // kT
			vec![h], // hists
			false // cyclic
		)
	}
	
	#[test]
	fn calc_bias() {
		let ds = build_hist_set(); // k = 10

		// 7th element -> x=7.5, x0=7.5
		assert_delta!(1.0, ds.calc_bias(7, 0), 0.00000001);

		// 8th element -> x=8.5, x0=7.5
		assert_delta!(0.13472233779, ds.calc_bias(8,0), 0.00000001);
		
		// 1st element -> x=0.5, x0=7.5. non-cyclic!
		assert_delta!(0.0, ds.calc_bias(0,0), 0.0000001);
	}

	#[test]
	fn calc_biascyclic() {
		let mut ds = build_hist_set();
		ds.cyclic = true;

		// 7th element -> x=7.5, x0=7.5
		assert_delta!(1.0, ds.calc_bias(7, 0), 0.00000001);

		// 8th element -> x=8.5, x0=7.5
		assert_delta!(0.13472233779, ds.calc_bias(8, 0), 0.00000001);
		

		// 1th element -> x=0.5, x0=7.5
		// cyclic flag makes bin 0 neighboring bin 9, so the distance is actually 2
		assert_delta!(0.00032942643, ds.calc_bias(0, 0), 0.00000001);

		// 2nd element -> x=1.5, x0=7.5
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
}