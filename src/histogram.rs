use std::fmt;

// One histogram
#[derive(Debug)]
pub struct Histogram {
	// offset of this histogram bins from the global histogram
	pub first: usize,

	// offset of the last element of the histogram. TODO required?
	pub last: usize,

	// total number of data points stored in the histogram
	pub num_points: u32,

	// histogram bins
	pub bins: Vec<f32>
}

impl Histogram {
	pub fn new(first: usize, last: usize, num_points: u32, bins: Vec<f32>) -> Histogram {
		Histogram {first, last, num_points, bins}
	}

	// Returns the value of a bin if the bin is present in this
	// histogram
	pub fn get_bin_count(&self, bin: usize) -> Option<f32> {
		if bin < self.first || bin > self.last {
			None
		} else {
			Some(self.bins[bin-self.first])
		}
	}
}

// a set of histograms

#[derive(Debug)]
pub struct HistogramSet {
	// number of histogram windows (number of simulations)
	pub num_windows: usize,

	// number of global histogram bins
	pub num_bins: usize,

	// min value of the histogram
	pub hist_min: f32,

	// max value of the histogram
	pub hist_max: f32,

	// width of a bin in unit of x
	pub bin_width: f32,

	// locations of biases
	pub bias_x0: Vec<f32>,

	// force constants of biases
	pub bias_fc: Vec<f32>,

	// value of kT
	pub kT: f32,

	// histogram for each window
	pub histograms: Vec<Histogram>
}

impl HistogramSet {
	
	pub fn new(num_bins: usize, bin_width: f32, hist_min: f32, hist_max: f32, bias_x0: Vec<f32>, bias_fc: Vec<f32>, kT: f32, histograms: Vec<Histogram>) -> HistogramSet {
		let num_windows = histograms.len();
		HistogramSet{num_windows, num_bins, bin_width, hist_min, hist_max, bias_x0, bias_fc, kT, histograms}
	}
}

impl fmt::Display for HistogramSet {
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

	fn build_hist() -> Histogram {
		Histogram{
			first: 5,
			last: 7,
			num_points: 5,
			bins: vec![1.0,1.0,3.0]
		}
	}

	#[test]
	fn get_bin_count() {
		let h = build_hist();
		assert_eq!(3.0, h.get_bin_count(7).unwrap());
		assert_eq!(3.0, h.get_bin_count(7).unwrap()); // twice for borrow
		assert_eq!(1.0, h.get_bin_count(5).unwrap());
		assert_eq!(None, h.get_bin_count(4));
		assert_eq!(None, h.get_bin_count(8));
	}
}