use rgsl::statistics;

// calculates the statistical inefficiency g of the given timeseries
// the quantity g can be thought of: N/g is the number of uncorrelated
// configurations in the timeseries, where samples are separated by
// the a multiple of g 
// For details, see "Chodera et al. (2007). Use of a Weighted Histogram Analysis
// Method for the Analysis of Simulated and Parallel Tempering Simulations, JCTC"
fn statistical_ineff(timeseries: &[f64]) -> f64 {
    let n = timeseries.len();
    let mean = statistics::mean(timeseries, 1, timeseries.len());
    let d_mean = timeseries.iter().map(|x| x-mean).collect::<Vec<f64>>();
    let cov = statistics::covariance(timeseries, 1, timeseries, 1, n);

    let mut g = 1.0;
    for t in 1..(n-1) {

        // normalized autocorr C(t) = (<x_n*x_(n+t)> - <x_n>^2) / (<x_n^2> - <x_n>^2)
        let tmp = d_mean[0..n-t].iter().zip(d_mean[t..n].iter()).map(|(x, y)| x*y);
        let c: f64 = tmp.map(|x| x+x).sum::<f64>() / (2.0 * (n as f64-t as f64)*cov);

        if c <= 0.0 {  // terminate at first 0 (autocorr gets noisy from here)
            break;
        }
        g = g + (2.0*c*(1.0-t as f64/n as f64))
    }
    if g < 1.0 {
        1.0
    } else {
        g
    }
}

// The autocorrelation time of a timeseries can be deduced from the
// `statistical_ineff` by (g-1)/2.0 
fn autocorrelation_time(g: f64) -> f64 {
    (g - 1.0) / 2.0
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead, BufReader};
    use std::fs::File;
    
    fn read_timeseries(filename: &str) -> Vec<f64> {
        let mut timeseries: Vec<f64> = Vec::new();
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(&file);
        for line in reader.lines() {
            let val = line.unwrap().split_whitespace()
                .collect::<Vec<&str>>()[1].parse::<f64>().unwrap();
            timeseries.push(val);
        }
        timeseries

    }

    #[test]
    fn statistical_ineff() {
        let timeseries = read_timeseries("example/1d_cyclic/COLVAR-2.5.xvg");
        let g = super::statistical_ineff(&timeseries);
        println!("{:?}", g);
        assert!((g - 3.859).abs() < 0.001)
    }

    #[test]
    fn autocorrelation_time() {
        let timeseries = read_timeseries("example/1d_cyclic/COLVAR-2.5.xvg");
        let g = super::statistical_ineff(&timeseries);
        let tau = super::autocorrelation_time(g);
        println!("{:?}", tau);
        assert!((tau - 1.430).abs() < 0.001)
    }

}