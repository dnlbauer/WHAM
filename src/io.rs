use super::histogram::Dataset;
use super::histogram::Histogram;
use super::Config;
use super::correlation_analysis::{statistical_ineff, autocorrelation_time};
use std::fs::OpenOptions;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader,BufWriter};
use k_B;
use std::path::Path;
use super::errors::*;
use f64;

// Returns the path to path2 relative to path1
// path1: "path/to/file.dat"
// path2: "another_file.dat"
// => result = path/to/another_file.dat
fn get_relative_path(path1: &str, path2: &str) -> String {
    let path1 = Path::new(path1);
    path1.parent().unwrap().join(path2).to_str().unwrap().to_string()
}

pub fn vprintln(s: String, verbose: bool) {
    if verbose {
        println!("{}", s);
    }
}

// Read input data into a histogram set by iterating over input files
// given in the metadata file. This generates at least one Dataset,
// or multiple Datasets if convdt is set in the config
pub fn read_data(cfg: &Config) -> Result<Vec<Dataset>> {
    let mut bias_pos: Vec<f64> = Vec::new();
    let mut bias_fc: Vec<f64> = Vec::new();
    let mut timeseries_lengths: Vec<usize> = Vec::new();
    let mut paths = Vec::new();

    // Boundaries of individual histograms if convdt is set.
    let dataset_boundaries: Vec<(f64, f64)> = get_convdt_boundaries(cfg.start, cfg.end, cfg.convdt);
    let num_datasets = dataset_boundaries.len();
    
    // for each timeseries, histograms are build for slices according to
    // start..convdt, start..2*convdt, ...
    let mut histograms =  vec![Vec::new(); dataset_boundaries.len()];

    let kT = cfg.temperature * k_B;
    let bin_width: Vec<f64> = (0..cfg.dimens).map(|idx| {
            (cfg.hist_max[idx] - cfg.hist_min[idx])/(cfg.num_bins[idx] as f64)
        }).collect();
    let num_bins: usize = cfg.num_bins.iter().product();
    let dimens_length = cfg.num_bins.clone();

    let f = File::open(&cfg.metadata_file).chain_err(|| "Failed to open metadata file")?;
    let buf = BufReader::new(&f);

    // read each metadata file line and parse it
    for (line_num,l) in buf.lines().enumerate() {
        let line = l.chain_err(|| "Failed to read line")?;

        // skip comments and empty lines
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let split: Vec<&str> = line.split_whitespace().collect();
        if split.len() < 1 + cfg.dimens * 2 {
            bail!(format!("Wrong number of columns in line {} of metadata file. Empty Line?", line_num+1));
        }

        // parse bias force constants and positions
        for val in split.iter().skip(1).take(cfg.dimens) {
            let pos = val.parse()
                .chain_err(|| format!("Failed to read bias position in line {} of metadata file", line_num+1))?;
            bias_pos.push(pos);
        }
        for val in split.iter().skip(1+cfg.dimens).take(cfg.dimens) {
            let fc = val.parse()
                .chain_err(|| format!("Failed to read bias fc in line {} of metadata file", line_num+1))?;
            bias_fc.push(fc);
        }

        // parse histogram data
        let path = get_relative_path(&cfg.metadata_file, split[0]);
        paths.push(path.clone());
        let (timeseries, timeseries_initial_lengths) = read_window_file(&path, cfg)
            .chain_err(|| format!("Failed to read time series from {}", &path))?;
        timeseries_lengths.push(timeseries_initial_lengths);

        for (idx, interval) in dataset_boundaries.iter().enumerate() {
            // build histogram for slice start.._stop
            let (start, stop) = interval;
            let timeseries_mask: Vec<bool> = (0..timeseries[0].len()).map(|i| {
                is_in_time_boundaries(timeseries[0][i], *start, *stop)
            }).collect();
            let hist = build_histogram_from_timeseries(&timeseries, &timeseries_mask, cfg);
            histograms[idx].push(hist);

            if (cfg.convdt == 0.00) || idx+1 == num_datasets {
                vprintln(format!("{}, {} data points added.",
                    &path, histograms[idx].last().unwrap().num_points), cfg.verbose);
                break
            }
        }
    }

    // Datasets are created from histograms.
    // Empty histograms result in an error when its the final dataset, and a warning otherwise.
    vprintln(format!("Generating {} datasets from histograms.", num_datasets), cfg.verbose);
    let datasets: Vec<Dataset> = histograms.into_iter().enumerate().map(|(dataset_idx, dataset_histograms)| {
        for (hs, path) in dataset_histograms.iter().zip(&paths) {
            if hs.num_points == 0 {
                let warning = format!("No data points for interval {}-{} in histogram boundaries: {}.",
                    dataset_boundaries[dataset_idx].0, dataset_boundaries[dataset_idx].1 ,&path);

                if dataset_idx+1 == num_datasets {
                    let warning = warning + " This is the final dataset.";
                    if cfg.ignore_empty {
                        eprintln!("{}", warning);
                    } else {
                        bail!(warning);
                    }
                } else {
                    eprintln!("{}", warning);
                }
            }
        }

        Ok(Dataset::new(num_bins, dimens_length.clone(), bin_width.clone(),
            cfg.hist_min.clone(), cfg.hist_max.clone(), bias_pos.clone(),
            bias_fc.clone(), kT, dataset_histograms, cfg.cyclic))
    }).collect::<Result<Vec<Dataset>>>().chain_err(|| "Failed to create datasets.")?;

    if datasets.is_empty() {
        bail!("No datasets created.")
    } else if datasets[0].histograms.is_empty() {
        bail!("Dataset has no associated data points.")
    } else {
        if datasets.len() > 1 {
            println!("Datasets:");
            println!("Dataset\t\tTime interval\t\tWindows\t\tN_total");
            for (idx, dataset) in datasets.iter().enumerate() {
                let n: u32 = dataset.histograms.iter().map(|h| h.num_points).sum();
                let mut stop = cfg.start+cfg.convdt*(idx+1) as f64;
                if stop > cfg.end {
                    stop = cfg.end;
                }
                println!("{:?}\t\t{:?}-{:?}\t\t{:?}\t\t{:?}", idx+1, cfg.start, stop, dataset.histograms.len(), n);
            }
        }

        let histograms = &datasets.last().unwrap().histograms;
        if cfg.uncorr {
            println!("Timeseries Correlation:");
            println!("Window\t\tN\t\tN_uncorr\tN/N_uncorr");
            for (idx, (n, h)) in timeseries_lengths.iter().zip(histograms.iter()).enumerate() {
                println!("{:?}\t\t{:?}\t\t{:?}\t\t{:.2}",
                    idx+1, n, h.num_points, h.num_points as f64 / *n as f64);
            }
            let total_n = timeseries_lengths.iter().sum::<usize>() as f64;
            let total_h = histograms.iter().map(|h| h.num_points).sum::<u32>() as f64;
            println!("\t\t\t\t\tTotal:\t{:.2}", total_h/total_n);
        }

        Ok(datasets)
    }
}

// builds a time boundaries for datasets from convdt, start and end
fn get_convdt_boundaries(start: f64, end: f64, convdt: f64) -> Vec<(f64, f64)> {
    if convdt == 0.0 {
        vec![(start, end)]
    } else {
        let intervals: usize = ((end - start) / convdt).ceil() as usize;
        (1..intervals+1).map(|i| {
            let interval_end = i as f64 * convdt + start;
            if interval_end > end {
                end
            } else {
                interval_end
            }
        }).map(|interval_end| { (start, interval_end) }).collect()
    }
}

// build a histogram from a timeseries
// mask is used to filter the timeseries for selected frames
fn build_histogram_from_timeseries(timeseries: &[Vec<f64>], mask: &[bool],
    cfg: &Config) -> Histogram {

    // total number of bins is the product of all dimensions length
    let total_bins = cfg.num_bins.iter().product();

    // bin width for each dimension: (max-min)/bins
    let bin_width: Vec<f64> = (0..cfg.dimens).map(|idx| {
        (cfg.hist_max[idx] - cfg.hist_min[idx])/(cfg.num_bins[idx] as f64)
    }).collect();

    // build histogram for slice start..convdt_stop
    let mut hist = vec![0.0; total_bins];
    for i in (0..timeseries[0].len()).filter(|i| mask[*i]) {
        let mut values: Vec<f64> = vec![f64::NAN; cfg.dimens+1];
        for j in 0..values.len() {
            values[j] = timeseries[j][i];
        }

        if is_in_hist_boundaries(&values[1..], cfg) {
            let bin_indeces: Vec<usize> = (0..cfg.dimens).map(|dimen: usize| {
                let val = values[dimen+1];
                ((val - cfg.hist_min[dimen]) / bin_width[dimen]) as usize
            }).collect();
            let index = flat_index(&bin_indeces, &cfg.num_bins);
            hist[index] += 1.0;
        }
    }

    let num_points: f64 = hist.iter().sum();
    Histogram::new(num_points as u32, hist)    
}

// transforms a multidimensional index into a one dimensional index
// indeces: multidimensional indeces
// lengths: length of the matrix in each dimension
// returns an index if the matrix is flattened to a one dimensional vector
// example for 3 dimensions N,M,O: idx = i_O + l_O*l_M*i_M + l_O*l_M*l_N*i_N
fn flat_index(indeces: &[usize], lengths: &[usize]) -> usize {
    indeces.iter().enumerate().map(|(i, idx)| {
        idx * lengths.iter().take(i).product::<usize>()
    }).sum()
}

// returns true if the values are inside the histogram boundaries defined by cfg
fn is_in_hist_boundaries(values: &[f64], cfg: &Config) -> bool {
    for dimen in 0..cfg.dimens {
        if values[dimen] < cfg.hist_min[dimen] || values[dimen] >= cfg.hist_max[dimen] {
            return false
        }
    }
    true
}

// returns true given time in inside the time boundaries defined by cfg
fn is_in_time_boundaries(time: f64, start: f64, end: f64) -> bool {
    if start <= time && time <= end {
        return true
    }
    false
}

// parse a time series file
fn read_window_file(window_file: &str, cfg: &Config) -> Result<(Vec<Vec<f64>>, usize)> {
    let mut timeseries: Vec<Vec<f64>> = read_timeseries(window_file, cfg)?;

    // filter the timeseries based on start/end parameters
    let time_series_mask: Vec<bool> = timeseries[0].iter()
        .map(|t| is_in_time_boundaries(*t, cfg.start, cfg.end)).collect();
    timeseries = timeseries.into_iter().map(|ts| {
        ts.into_iter().zip(time_series_mask.iter()).filter_map(|(val, mask)| {
            if *mask {
                Some(val)
            } else {
                None
            }
        }).collect()
    }).collect::<Vec<Vec<f64>>>();

    let timeseries_inital_length = timeseries[0].len();
    if cfg.uncorr {
        timeseries = uncorrelate(timeseries, cfg);
    }

    if timeseries[0].is_empty() && !cfg.ignore_empty {
        bail!("Time series is empty")
    }

    Ok((timeseries, timeseries_inital_length))
}

// Read a multidimensional timeseries
// The resulting vector contains one vector per dimension
fn read_timeseries(window_file: &str, cfg: &Config) -> Result<Vec<Vec<f64>>> {
    let f = File::open(window_file)
        .chain_err(|| format!("Failed to open sample data file {}.", window_file))?;
    let mut buf = BufReader::new(&f);

    let mut timeseries = vec![Vec::new(); cfg.dimens+1];

    // read and parse each timeseries line
    let mut line = String::new();
    let mut linecount = 0;
    while buf.read_line(&mut line).chain_err(|| "Failed to read line")? > 0 {
        linecount += 1;

        // skip comments and empty lines
        if line.starts_with('#') || line.starts_with('@') || line.is_empty() {
            line.clear();
            continue;
        }

        {
            let split: Vec<&str> = line.split_whitespace().collect();
            if split.len() < cfg.dimens+1 {
                bail!(format!("Wrong number of columns in line {} of window file {}. Empty Line?.", linecount, window_file));
            }

            for i in 0..cfg.dimens+1 {
                timeseries[i].push(split[i].parse::<f64>()
                    .chain_err(|| format!("Failed to parse line {} of window file {}.", linecount, window_file))?

                );
            }
        }
        
        line.clear();
    }
    Ok(timeseries)
}


// calculates the inefficiency for every collective variable
// filters the timeseries based on the highest inefficiency
fn uncorrelate(timeseries: Vec<Vec<f64>>, cfg: &Config) -> Vec<Vec<f64>> {
    // calculate inefficiencies and find the highest one
    let gs: Vec<f64> = timeseries[1..].iter().map(|ts| statistical_ineff(ts)).collect();
    let mut max_g = 1.0;
    for g in gs {
        if g > max_g {
            max_g = g;
        }
    }

    // round g up
    let mut trunc_g = max_g.trunc() as usize;
    if (trunc_g as f64 - max_g).abs() > 0.000_000_000_1 {
        trunc_g += 1;
    }

    // filter correlated samples from timeseries
    let prev_len = timeseries[0].len();
    let timeseries = timeseries.into_iter().map(|ts| {
        ts.into_iter().step_by(trunc_g).collect::<Vec<f64>>()
    }).collect::<Vec<Vec<f64>>>();

    let new_len = timeseries[0].len();
    if cfg.verbose {
        let tau = autocorrelation_time(max_g)* (timeseries[0][1]-timeseries[0][0]);
        vprintln(format!("{:?}/{:?} samples are uncorrelated. {:?} samples removed from timeseries (tau={:.5})", new_len, prev_len, prev_len-new_len, tau), true);
    }
    timeseries
}

// Write WHAM calculation results to out_file.
pub fn write_results(out_file: &str, append: bool, ds: &Dataset, free: &[f64],
    free_std: &[f64], prob: &[f64], prob_std: &[f64], index: Option<usize>) -> Result<()> {

    if !append && Path::new(out_file).exists() {
        std::fs::remove_file(out_file).chain_err(|| "Failed to delete file.")?;
    }
    let output = OpenOptions::new().write(true)
        .append(true)
        .create(true)
        .open(out_file)
        .chain_err(|| format!("Failed to create file with path {}", out_file))?;
    let mut buf = BufWriter::new(output);

    let header: String = (0..ds.dimens_lengths.len()).map(|d| format!("coord{}", d+1))
        .collect::<Vec<String>>().join("    ");
    if let Some(index) = index {
        writeln!(buf, "#Dataset {}", index).unwrap();
    }
    writeln!(buf, "#{}    Free Energy    +/-    Probability    +/-", header).unwrap();

    for bin in 0..free.len() {
        let coords = ds.get_coords_for_bin(bin);
        let coords_str: String = coords.iter().map(|c| {format!("{:8.6}    ", c)})
            .collect::<Vec<String>>().join("\t");
        writeln!(buf, "{}{:8.6}    {:8.6}    {:8.6}    {:8.6}", coords_str,
            free[bin], free_std[bin], prob[bin], prob_std[bin])
            .chain_err(|| "Failed to write to file.")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    fn cfg() -> Config {
        Config {
            metadata_file: "example/1d_cyclic/metadata.dat".to_string(),

            hist_min: vec![-3.14],
            hist_max: vec![3.14],
            num_bins: vec![10],
            dimens: 1,
            verbose: false,
            tolerance: 0.0,
            max_iterations: 0,
            temperature: 300.0,
            cyclic: false,
            output: "qwert".to_string(),
            bootstrap: 0,
            bootstrap_seed: 1234,
            start: 0.0,
            end: 1e+20,
            uncorr: false,
            convdt: 0.0,
            ignore_empty: false,
        }
    }

    #[test]
    fn read_window_file() {
        let f = "example/1d_cyclic/COLVAR+0.0.xvg";
        let cfg = cfg();
        let (timeseries, timeseries_inital_length) = super::read_window_file(&f, &cfg).unwrap();
        let mask = vec![true; timeseries[0].len()];
        let h = build_histogram_from_timeseries(&timeseries, &mask, &cfg);
        println!("{:?}", h);
        assert_eq!(5000, timeseries_inital_length);
        assert_eq!(5000, h.num_points);
        assert_approx_eq!(0.0, h.bins[2]);
        assert_approx_eq!(11.0, h.bins[3]);
        assert_approx_eq!(2236.0, h.bins[4]);
        assert_approx_eq!(2714.0, h.bins[5]);
        assert_approx_eq!(39.0, h.bins[6]);
        assert_approx_eq!(0.0, h.bins[7]);
    }

    #[test]
    fn read_timeseries() {
        let f = "example/1d_cyclic/COLVAR+0.0.xvg";
        let cfg = cfg();
        let ts = super::read_timeseries(&f, &cfg).unwrap();
        let expected = [
            -0.153_145,
            -0.377_860,
            0.010_992,
            0.123_074,
            0.108_291,
            0.261_607,
        ];
        assert!(ts.len() == 2);
        assert!(ts[0].len() == 5000);
        println!("{:?}", ts);
        for (actual, expected) in ts[1].iter().zip(expected.iter()) {
            assert!((actual-expected).abs() < 0.001, format!("{:?} != {:?}", actual, expected));
        }
    }

    #[test]
    fn read_data() {
        let cfg = cfg();
        let ds = &super::read_data(&cfg).unwrap()[0];
        println!("{:?}", ds);
        assert_eq!(25, ds.num_windows);
        assert_eq!(cfg.num_bins.len(), ds.dimens_lengths.len());
        assert_eq!(cfg.num_bins[0], ds.dimens_lengths[0]);
        assert_approx_eq!(cfg.temperature * k_B, ds.kT);
        assert_eq!(25, ds.histograms.len())
    }

    #[test]
    fn read_data_empty() {
        let mut cfg = cfg();
        cfg.metadata_file = "tests/data/metadata_convdt.dat".to_string();
        cfg.start = 2.5;
        cfg.end = 9.0;
        cfg.ignore_empty = false;

        // should throw an error since one first timeseries ends at 2
        let ds = super::read_data(&cfg);
        if ds.is_ok() {
            panic!()
        }

        // should not throw an error because ignore_empty is set
        cfg.ignore_empty = true;
        let ds = super::read_data(&cfg);
        if ds.is_err() {
            panic!()
        }
    }

    // test if convdt results in correct parsing
    // 6 timeseries are loaded ranging from:
    // 1. 0-10, 500 datapoints
    // 2. 0-2,  100 datapoints
    // 3. 0-5,  250 datapoints
    // 4. 5-10, 250 datapoints
    // 5. 7-10, 150 datapoints
    // 6  2-7,  250 datapoints
    #[test]
    fn read_data_convdt() {
        let mut cfg = cfg();
        cfg.metadata_file = "tests/data/metadata_convdt.dat".to_string();
        cfg.convdt = 2.0;
        cfg.start = 0.0;
        cfg.end = 9.0;
        let dss = super::read_data(&cfg).unwrap();
        assert_eq!(5, dss.len());

        for ds in &dss {
            assert_eq!(6, ds.num_windows);
            assert_eq!(6, ds.histograms.len());
        }

        let hist_points: Vec<u32> = dss.iter().map(|ds| {
            ds.histograms.iter().map(|h| h.num_points).sum()
        }).collect();

        let expected_hist_points = vec![
            300,  // 0-2: 100+100+100+0+0
            600,  // 0-4: 200+100+200+0+0+100
            900,  // 0-6: 300+100+250+50+0+200
            1200, // 0-8: 400+100+250+150+50+250
            1350, // 0-9: 450+100+250+200+100+250
        ];
        for (expected, actual) in expected_hist_points.iter().zip(hist_points.iter()) {
            assert_eq!(expected, actual);
        }
    }

    // test convdt with a single time series
    #[test]
    fn read_data_convdt_single() {
        let mut cfg = cfg();
        cfg.metadata_file = "tests/data/metadata_convdt_single.dat".to_string();
        cfg.convdt = 2.0;
        cfg.start = 0.0;
        cfg.end = 9.0;
        let dss = super::read_data(&cfg).unwrap();
        assert_eq!(5, dss.len());

        for ds in &dss {
            assert_eq!(1, ds.num_windows);
            assert_eq!(1, ds.histograms.len());
        }

        let hist_points: Vec<u32> = dss.iter().map(|ds| {
            ds.histograms.iter().map(|h| h.num_points).sum()
        }).collect();

        let expected_hist_points = vec![
            0,   // 0-2
            100, // 0-4
            200, // 0-6
            250, // 0-8
            250, // 0-9
        ];
        for (expected, actual) in expected_hist_points.iter().zip(hist_points.iter()) {
            assert_eq!(expected, actual);
        }
    }

    #[test]
    fn get_relative_path() {
        let path1 = "path/to/some_file.dat";
        let path2 = "another_file.dat";
        let path3 = "subfolder/another_file.dat";
        let relative2 = super::get_relative_path(&path1, &path2);
        assert_eq!("path/to/another_file.dat" ,relative2);
        let relative3 = super::get_relative_path(&path1, &path3);
        assert_eq!("path/to/subfolder/another_file.dat" ,relative3);
    }

    #[test]
    fn is_in_time_boundaries() {
        let start = 10.0;
        let end = 20.0;
        assert!(super::is_in_time_boundaries(15.0, start, end));
        assert!(super::is_in_time_boundaries(10.0, start, end));
        assert!(super::is_in_time_boundaries(20.0, start, end));
        assert!(!super::is_in_time_boundaries(9.9999999, start, end));
        assert!(!super::is_in_time_boundaries(20.000001, start, end));   
    }

    #[test]
    fn get_convdt_boundaries() {
        let test = super::get_convdt_boundaries(10.0, 20.0, 10.0);
        println!("{:?}", test);
        assert!(test.len() == 1);
        assert_approx_eq!(test[0].0, 10.0);
        assert_approx_eq!(test[0].1, 20.0);

        let test = super::get_convdt_boundaries(10.0, 20.0, 5.0);
        println!("{:?}", test);
        assert!(test.len() == 2);
        assert_approx_eq!(test[0].0, 10.0);
        assert_approx_eq!(test[0].1, 15.0);
        assert_approx_eq!(test[1].0, 10.0);
        assert_approx_eq!(test[1].1, 20.0);

        let test = super::get_convdt_boundaries(5.0, 30.0, 10.0);
        println!("{:?}", test);
        assert!(test.len() == 3);
        assert_approx_eq!(test[0].0, 5.0);
        assert_approx_eq!(test[0].1, 15.0);
        assert_approx_eq!(test[1].0, 5.0);
        assert_approx_eq!(test[1].1, 25.0);
        assert_approx_eq!(test[2].0, 5.0);
        assert_approx_eq!(test[2].1, 30.0);
    }
}