use std::process::Command;
use std::fs;

#[cfg(test)]
mod integration {
    use Command;
    use fs;

    #[test]
    fn no_args() {
        let output = Command::new("./target/debug/wham")
            .output()
            .expect("failed to execute process");
        let expected_error = "The following required arguments were not provided";
        let out = String::from_utf8_lossy(&output.stderr);
        assert!(out.contains(expected_error));
    }

    #[test]
    fn unparseable_bias_pos() {
        let output = Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--min", "-3.0", "--max", "3.0", "-T", "300"])
            .args(&["-f", "tests/data/metadata_unparseable1.dat"])
            .args(&["-o", "/dev/null"])
            .output()
            .expect("failed to execute process");
        let output = String::from_utf8_lossy(&output.stderr);
        println!("{}", output);
        assert!(output.to_string().contains(
            "Failed to read bias position in line 1"
        ));
    }

    #[test]
    fn unparseable_bias_fc() {
        let output = Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--min", "-3.0", "--max", "3.0", "-T", "300"])
            .args(&["-f", "tests/data/metadata_unparseable2.dat"])
            .args(&["-o", "/dev/null"])
            .output()
            .expect("failed to execute process");
        let output = String::from_utf8_lossy(&output.stderr);
        println!("{}", output);
        assert!(output.to_string().contains(
            "Failed to read bias fc in line 1"
        ));
    }

    #[test]
    fn no_convergence() {
        let output = Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--min", "-3.0", "--max", "3.0", "-T", "300"])
            .args(&["--iterations", "10"])
            .args(&["-f", "example/1d_cyclic/metadata.dat"])
            .args(&["-o", "/dev/null"])
            .output()
            .expect("failed to execute process");
        let output = String::from_utf8_lossy(&output.stderr);
        println!("{}", output);
        assert!(output.to_string().contains(
            "Error: WHAM not converged! (max iterations reached)"
        ));
    }

    #[test]
    fn unparseable_timeseries() {
        let output = Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--min", "-3.0", "--max", "3.0", "-T", "300"])
            .args(&["-f", "tests/data/metadata_unparseable3.dat"])
            .args(&["-o", "/dev/null"])
            .output()
            .expect("failed to execute process");
        let output = String::from_utf8_lossy(&output.stderr);
        println!("{}", output);
        assert!(output.to_string().contains(
            "Failed to parse line"
        ));
    }

    #[test]
    fn unparseable_timeseries_empty() {
        let output = Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--min", "-3.0", "--max", "3.0", "-T", "300"])
            .args(&["-f", "tests/data/metadata_unparseable4.dat"])
            .args(&["-o", "/dev/null"])
            .output()
            .expect("failed to execute process");
        let output = String::from_utf8_lossy(&output.stderr);
        println!("{}", output);
        assert!(output.to_string().contains(
            "Wrong number of columns in line"
        ));
    }

    #[test]
    fn skip_rows() {
        let output = Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--max", "3.14", "--min", "-3.14", "-T", "300", "--cyclic"])
            .args(&["-f", "example/1d_cyclic/metadata.dat"])
            .args(&["-o", "/tmp/wham_test_1d_cyclic.out"])
            .args(&["--start", "50", "--end", "60"])
            .args(&["-v"])
            .output()
            .expect("failed to execute process");

        let output = String::from_utf8_lossy(&output.stdout);
        println!("{}", output);
        assert!(output.to_string().contains(
            "25 windows, 12489 datapoints"
        ));
    }

}