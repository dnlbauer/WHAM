use super::histogram::Dataset;
use super::histogram::Histogram;
use super::Config;
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
// given in the metadata file
pub fn read_data(cfg: &Config) -> Result<Dataset> {
	let mut bias_pos: Vec<f64> = Vec::new();
	let mut bias_fc: Vec<f64> = Vec::new();
    let mut histograms: Vec<Histogram> = Vec::new();

	let kT = cfg.temperature * k_B;
    let bin_width: Vec<f64> = (0..cfg.dimens).map(|idx| {
            (cfg.hist_max[idx] - cfg.hist_min[idx])/(cfg.num_bins[idx] as f64)
        }).collect();
    let num_bins = cfg.num_bins.iter().fold(1, |state, &bins| state*bins);
    let dimens_length = cfg.num_bins.clone();

    let f = File::open(&cfg.metadata_file).chain_err(|| "Failed to open metadata file")?;
    let buf = BufReader::new(&f);

    // read each metadata file line and parse it
    for (line_num,l) in buf.lines().enumerate() {
    	let line = l.chain_err(|| "Failed to read line")?;

        // skip comments and empty lines
        if line.starts_with("#") || line.len() == 0 {
    		continue;
    	}

    	let split: Vec<&str> = line.split_whitespace().collect();
        if split.len() < 1 + cfg.dimens * 2 {
            bail!(format!("Wrong number of columns in line {} of metadata file. Empty Line?", line_num+1));
        }

        // parse histogram data
        let path = get_relative_path(&cfg.metadata_file, split[0]);
        let h = read_window_file(&path, cfg)
            .chain_err(|| format!("Failed to parse process data file {}", &path))?;
        if h.num_points == 0 {
            bail!(format!("No data points in histogram boundaries: {}", &path))
        }
        histograms.push(h);
        vprintln(format!("{}, {} data points added.", &path,
                         histograms.last().unwrap().num_points), cfg.verbose);

        // parse bias force constants and positions
        for i in 1..cfg.dimens+1 {
            let pos = split[i].parse()
                .chain_err(|| format!("Failed to read bias position in line {} of metadata file", line_num+1))?;
            bias_pos.push(pos);
        }
        for i in (1+cfg.dimens)..(1+2*cfg.dimens) {
            let fc = split[i].parse()
                .chain_err(|| format!("Failed to read bias fc in line {} of metadata file", line_num+1))?;
            bias_fc.push(fc);
        }
    }
    
    if histograms.len() > 0 {
        Ok(Dataset::new(num_bins, dimens_length, bin_width, cfg.hist_min.clone(), cfg.hist_max.clone(), bias_pos, bias_fc, kT, histograms, cfg.cyclic))
    } else {
        bail!("Histogram has no datapoints.")
    }
}

// transforms a multidimensional index into a one dimensional index
// indeces: multidimensional indeces
// lengths: length of the matrix in each dimension
// returns an index if the matrix is flattened to a one dimensional vector
// example for 3 dimensions N,M,O: idx = i_O + l_O*l_M*i_M + l_O*l_M*l_N*i_N
fn flat_index(indeces: &Vec<usize>, lengths: &Vec<usize>) -> usize {
    let mut idx = 0;
    for i in 0..indeces.len() {
        idx += indeces[i]*lengths[0..i].iter()
            .fold(1, |state, &l| { state * l });
    }
    idx
}

// returns true if the values are inside the histogram boundaries defined by cfg
fn is_in_hist_boundaries(values: &[f64], cfg: &Config) -> bool {
    for dimen in 0..cfg.dimens {
        if values[dimen] < cfg.hist_min[dimen] || values[dimen] > cfg.hist_max[dimen] {
            return false
        }
    }
    true
}

// returns true given time in inside the time boundaries defined by cfg
fn is_in_time_boundaries(time: f64, cfg: &Config) -> bool {
    if cfg.start <= time && time <= cfg.end {
        return true
    }
    false
}

// parse a timeseries file into a histogram
fn read_window_file(window_file: &str, cfg: &Config) -> Result<Histogram> {
	let f = File::open(window_file)
        .chain_err(|| format!("Failed to open sample data file {}.", window_file))?;
    let mut buf = BufReader::new(&f);

    // total number of bins is the product of all dimensions length
    let total_bins = cfg.num_bins.iter().fold(1, |s, &x| { s*x });
    let mut hist = vec![0.0; total_bins];

    // bin width for each dimension: (max-min)/bins
    let bin_width: Vec<f64> = (0..cfg.dimens).map(|idx| {
            (cfg.hist_max[idx] - cfg.hist_min[idx])/(cfg.num_bins[idx] as f64)
        }).collect();

    // read and parse each timeseries line
    let mut line = String::new();
    let mut linecount = 0;
    while buf.read_line(&mut line).chain_err(|| "Failed to read line")? > 0 {
        linecount += 1;
        // skip comments and empty lines
        if line.starts_with("#") || line.starts_with("@") || line.len() == 0 {
            line.clear();
            continue;
        }

        {

            let split: Vec<&str> = line.split_whitespace().collect();
            if split.len() < cfg.dimens+1 {
                bail!(format!("Wrong number of columns in line {} of window file {}. Empty Line?.", linecount, window_file));
            }

            let mut values: Vec<f64> = vec![f64::NAN; cfg.dimens+1];
            for i in 0..values.len() {
                values[i] = split[i].parse::<f64>()
                    .chain_err(|| format!("Failed to parse line {} of window file {}.", linecount, window_file))?;
            }

            if is_in_hist_boundaries(&values[1..], cfg) && is_in_time_boundaries(values[0], cfg) {
                let bin_indeces = (0..cfg.dimens).map(|dimen: usize| {
                    let val = values[dimen+1];
                    ((val - cfg.hist_min[dimen]) / bin_width[dimen]) as usize
                }).collect();
                let index = flat_index(&bin_indeces, &cfg.num_bins);
                hist[index] += 1.0;
            }
        }
        line.clear();
    }

    let num_points: f64 = hist.iter().sum();
    Ok(Histogram::new(num_points as u32, hist))
}

// Write WHAM calculation results to out_file.
pub fn write_results(out_file: &str, ds: &Dataset, free: &Vec<f64>, free_std: &Vec<f64>, prob: &Vec<f64>, prob_std: &Vec<f64>) -> Result<()> {
    let output = File::create(out_file)
        .chain_err(|| format!("Failed to create file with path {}", out_file))?;
    let mut buf = BufWriter::new(output);

    let header: String = (0..ds.dimens_lengths.len()).map(|d| format!("coord{}", d+1))
        .collect::<Vec<String>>().join("    ");
    writeln!(buf, "#{}    {}    {}    {}    {}", header, "Free Energy", "+/-", "Probability", "+/-");

    for bin in 0..free.len() {
        let coords = ds.get_coords_for_bin(bin);
        let coords_str: String = coords.iter().map(|c| {format!("{:8.6}    ", c)})
            .collect::<Vec<String>>().join("\t");
        writeln!(buf, "{}{:8.6}    {:8.6}    {:8.6}    {:8.6}", coords_str, free[bin], free_std[bin], prob[bin], prob_std[bin])
            .chain_err(|| "Failed to write to file.")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> Config {
        Config {
            metadata_file: "example/1d/metadata.dat".to_string(),
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
            start: 0.0,
            end: 1e+20
        }
    }

    #[test]
    fn read_window_file() {
        let f = "example/1d/COLVAR+0.0.xvg";
        let cfg = cfg();
        let h = super::read_window_file(&f, &cfg).unwrap();
        println!("{:?}", h);
        assert_eq!(5000, h.num_points);
        assert_eq!(0.0, h.bins[2]);
        assert_eq!(11.0, h.bins[3]);
        assert_eq!(2236.0, h.bins[4]);
        assert_eq!(2714.0, h.bins[5]);
        assert_eq!(39.0, h.bins[6]);
        assert_eq!(0.0, h.bins[7]);
    }


    #[test]
    fn read_data() {
        let cfg = cfg();
        let ds = super::read_data(&cfg).unwrap();
        println!("{:?}", ds);
        assert_eq!(25, ds.num_windows);
        assert_eq!(cfg.num_bins.len(), ds.dimens_lengths.len());
        assert_eq!(cfg.num_bins[0], ds.dimens_lengths[0]);
        assert_eq!(cfg.temperature * k_B, ds.kT);
        assert_eq!(25, ds.histograms.len())
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
}