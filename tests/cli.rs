mod command;

#[cfg(test)]
mod integration {
    use super::command::get_command;

    #[test]
    fn no_args() {
        let output = get_command()
            .output()
            .expect("failed to execute process");
        let expected_error = "The following required arguments were not provided";
        let out = String::from_utf8_lossy(&output.stderr);
        assert!(out.contains(expected_error));
    }

    #[test]
    fn unparseable_bias_pos() {
        let output = get_command()
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
        let output = get_command()
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
        let output = get_command()
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
        let output = get_command()
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
        let output = get_command()
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
        let output = get_command()
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

    #[test]
    fn convdt_needs_start_end() {
        let output = get_command()
            .args(&["--bins", "100", "--min", "-3.0", "--max", "3.0", "-T", "300"])
            .args(&["-f", "tests/data/metadata_unparseable1.dat"])
            .args(&["-o", "/dev/null"])
            .args(&["--convdt", "100"])
            .output()
            .expect("failed to execute process");
        let output = String::from_utf8_lossy(&output.stderr);
        println!("{}", output);
        assert!(output.to_string().contains(
            "--convdt requires --start and --end to be set"
        ));
    }

}