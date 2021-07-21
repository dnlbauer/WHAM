mod command;

#[cfg(test)]
mod integration {

    use std::process::Command;
    use std::fs;
    use super::command::get_command;
    use std::fs::OpenOptions;
    use std::io::prelude::*;

    #[test]
    fn wham_1d_cyclic() {
        let output_file = "/tmp/wham_test_1d_cyclic.out";
        get_command()
            .args(&["--bins", "100", "--max", "pi", "--min", "-pi", "-T", "300", "--cyclic"])
            .args(&["--bt", "100", "--seed", "1234"])
            .args(&["-f", "example/1d_cyclic/metadata.dat"])
            .args(&["-o", output_file])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata(output_file).is_ok());
        let output = Command::new("diff")
            .arg(output_file)
            .arg("example/1d_cyclic/wham.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
        std::fs::remove_file(output_file).unwrap();
    }

    #[test]
    fn wham_1d_cyclic_uncorrelated() {
        let output_file = "/tmp/wham_test_1d_cyclic.out";
        get_command()
            .args(&["--bins", "100", "--max", "pi", "--min", "-pi", "-T", "300", "--cyclic", "--uncorr"])
            .args(&["--seed", "1234"])
            .args(&["-f", "example/1d_cyclic/metadata.dat"])
            .args(&["-o", output_file])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata("/tmp/wham_test_1d_cyclic.out").is_ok());
        let output = Command::new("diff")
            .arg(output_file)
            .arg("example/1d_cyclic/wham_uncorrelated.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
        std::fs::remove_file(output_file).unwrap();
    }

    #[test]
    // Test if convdt runs have the same result as normal runs
    fn wham_convdt() {
        // run wham with convdt
        let output_file = "/tmp/wham_test_convdt.out";
        get_command()
            .args(&["--bins", "10", "--max", "pi", "--min", "-pi", "-T", "300", "--cyclic"])
            .args(&["--seed", "1234"])
            .args(&["--start", "0", "--end", "100"])
            .args(&["--convdt", "10"])
            .args(&["-f", "example/1d_cyclic/metadata.dat"])
            .args(&["-o", output_file])
            .output()
            .expect("failed to execute process");
        assert!(fs::metadata(output_file).is_ok());

        // run wham for individual sets
        for i in (10..110).step_by(10) {
            let output_file_single = format!("/tmp/wham_test_convdt_{}.out", i);
            get_command()
                .args(&["--bins", "10", "--max", "pi", "--min", "-pi", "-T", "300", "--cyclic"])
                .args(&["--seed", "1234"])
                .args(&["--start", "0", "--end", &i.to_string()])
                .args(&["-f", "example/1d_cyclic/metadata.dat"])
                .args(&["-o", &output_file_single])
                .output()
                .expect("failed to execute process");
            assert!(fs::metadata(output_file_single).is_ok());
        }

        // combine individual runs
        let output_combined = "/tmp/wham_test_convdt_combined.out";
        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .open(output_combined)
            .unwrap();
        for i in 1..11 {
            let output_file_single = format!("/tmp/wham_test_convdt_{}.out", i*10);
            println!("{}", output_file_single);
            file.write_all(format!("#Dataset {}\n", i-1).as_bytes()).unwrap();
            file.write_all(fs::read_to_string(output_file_single.clone()).unwrap().as_bytes()).unwrap();
            std::fs::remove_file(output_file_single).unwrap();
        }

        // compare combined runs with single run
        let output = Command::new("diff")
            .arg(output_file)
            .arg(output_combined)
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);

        std::fs::remove_file(output_combined).unwrap();
        std::fs::remove_file(output_file).unwrap();
        
    }

    #[test]
    #[ignore] // expensive
    fn wham_2d_cyclic() {
        let output_file = "/tmp/wham_test_2d_cyclic.out";
        get_command()
            .args(&["--bins", "100,100", "--max", "pi,pi", "--min", "-pi,-pi", "-T", "300", "--cyclic"])
            .args(&["-f", "example/2d_cyclic/metadata.dat"])
            .args(&["-o", output_file])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata(output_file).is_ok());
        let output = Command::new("diff")
            .arg(output_file)
            .arg("example/2d_cyclic/wham.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
        std::fs::remove_file(output_file).unwrap();
    }
}
