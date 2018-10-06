use super::histogram::HistogramSet;
use super::histogram::Histogram;
use super::Config;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use k_B;
use std::process;
use std::option::Option;
use std::path::Path;

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
pub fn read_data(cfg: &Config) -> Option<HistogramSet> {
	let mut bias_x0: Vec<f32> = Vec::new();
	let mut bias_fc: Vec<f32> = Vec::new();
    let mut histograms: Vec<Histogram> = Vec::new();
	let kT = cfg.temperature * k_B;
    let f = File::open(&cfg.metadata_file).unwrap_or_else(|x| {
        eprintln!("Failed to read metadata from {}. {}", &cfg.metadata_file, x);
        process::exit(1)
        });
    let buf = BufReader::new(&f);

    for l in buf.lines() {
    	let line = l.unwrap();
        // skip comments and empty lines
        if line.starts_with("#") || line.len() == 0 {
    		continue;
    	}
    	let (path, x0, k) = scan_fmt!(&line, "{} {} {}", String, f32, f32);
        
        let path = get_relative_path(&cfg.metadata_file, &path.unwrap());
        match read_window_file(&path, cfg) {
            Some(h) => {
                histograms.push(h);
                vprintln(format!("{}, {} data points added.", &path, histograms.last().unwrap().num_points), cfg.verbose);
            },
            None => {
                eprintln!("No data points inside histogram boundaries: {}", &path);
                continue; // goto next iteration and skip adding bias values
            }
        }

    	bias_x0.push(x0.unwrap_or_else(|| {
            eprintln!("Failed to read coordinate from: {}", &line);
            process::exit(1)
        }));
    	bias_fc.push(k.unwrap_or_else(|| {
            eprintln!("Failed to read bias value from: {}", &line);
            process::exit(1)
        }));

        
    }
    
    if histograms.len() > 0 {
        let bin_width = (cfg.hist_max - cfg.hist_min)/(cfg.num_bins as f32);
        Some(HistogramSet::new(cfg.num_bins, bin_width, cfg.hist_min, cfg.hist_max, bias_x0, bias_fc, kT, histograms)) 
    } else {
        None
    }
}

// parse a timeseries file into a histogram
fn read_window_file(window_file: &str, cfg: &Config) -> Option<Histogram> {
	let f = File::open(window_file).unwrap_or_else(|x| {
        eprintln!("Failed to read sample data from {}. {}", window_file, x);
        process::exit(1)
    });
    let buf = BufReader::new(&f);
    
    let mut global_hist = vec![0.0; cfg.num_bins];
    let bin_width = (cfg.hist_max - cfg.hist_min)/(cfg.num_bins as f32);
    for l in buf.lines() {
    	let line = l.unwrap();
        // skip comments and empty lines
        if line.starts_with("#") || line.starts_with("@") || line.len() == 0 {
    		continue;
    	}

    	let (_, x) = scan_fmt!(&line, "{} {}", f32, f32);

    	match x {
    		Some(x) => {
    			if x > cfg.hist_min && x < cfg.hist_max {
    				let bin_ndx = ((x-cfg.hist_min) / bin_width) as usize;
                    global_hist[bin_ndx] += 1.0;
    			} 
    		}
    		None => { 
                eprintln!("{}, Failed to read datapoint from line: {}", &window_file, &line);
                process::exit(1);
            }
    	}
    }

    let mut max_bin: usize = 0;
    let mut min_bin: usize =(cfg.num_bins-1) as usize;
    for bin in 0..global_hist.len() {
        if global_hist[bin] != 0.0 && bin > max_bin {
            max_bin = bin;
        }
        if global_hist[bin] != 0.0 && bin < min_bin {
            min_bin = bin;
        }
    }

    if (max_bin == min_bin && global_hist[max_bin] == 0.0) || max_bin < min_bin {
        None // zero length histogram
    } else {
        // trim global hist to save memory
        global_hist.truncate(max_bin+1);
        global_hist.drain(..min_bin);
        
        let num_points: f32 = global_hist.iter().sum();
        Some(Histogram::new(min_bin, max_bin, num_points as u32, global_hist))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> Config {
        Config{
            metadata_file: "tests/data/metadata.dat".to_string(),
            hist_min: 0.0,
            hist_max: 3.0,
            num_bins: 30,
            verbose: false,
            tolerance: 0.0,
            max_iterations: 0,
            temperature: 300.0
        }
    }

    #[test]
    fn read_window_file() {
        let f = "tests/data/window_0.0.dat";
        let cfg = cfg();
        let h = super::read_window_file(&f, &cfg).unwrap();
        println!("{:?}", h);
        assert_eq!(1, h.first);
        assert_eq!(6, h.last);
        assert_eq!(11, h.num_points);
        assert_eq!(2.0, h.bins[0]);
        assert_eq!(2.0, h.bins[2]);
        assert_eq!(2.0, h.bins[4]);
        assert_eq!(1.0, h.bins[5]);      
    }


    #[test]
    fn read_data() {
        let cfg = cfg();
        let hs = super::read_data(&cfg);
        assert!(hs.is_some());
        let hs = hs.unwrap();
        println!("{:?}", hs);
        assert_eq!(2, hs.num_windows);
        assert_eq!(cfg.num_bins, hs.num_bins);
        assert_eq!(cfg.hist_min, hs.hist_min);
        assert_eq!(cfg.hist_max, hs.hist_max);
        let expected_bin_width = (cfg.hist_max - cfg.hist_min)/cfg.num_bins as f32;
        assert_eq!(expected_bin_width, hs.bin_width);  
        assert_eq!(vec![0.0, 1.0], hs.bias_x0);  
        assert_eq!(vec![100.0, 200.0], hs.bias_fc);
        assert_eq!(cfg.temperature * k_B, hs.kT);
        assert_eq!(2, hs.histograms.len())
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