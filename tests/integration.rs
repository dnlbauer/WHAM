use std::process::Command;
use std::fs;

#[cfg(test)]
mod integration {
    use Command;
    use fs;

    #[test]
    fn calling_wham_without_args() {
        let output = Command::new("./target/debug/wham")
            .output()
            .expect("failed to execute process");
        let expected_error = "The following required arguments were not provided";
        let out = String::from_utf8_lossy(&output.stderr);
        assert!(out.contains(expected_error));
    }

    #[test]
    fn calling_wham_1d() {;
        Command::new("./target/debug/wham")
            .args(&["--bins", "100", "--max", "3.14", "--min", "-3.14", "-T", "300", "--cyclic"])
            .args(&["-f", "example/1d/metadata.dat"])
            .args(&["-o", "/tmp/wham_test_1d.out"])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata("/tmp/wham_test_1d.out").is_ok());
        let output = Command::new("diff")
            .arg("/tmp/wham_test_1d.out")
            .arg("example/1d/wham.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
    }

    #[test]
    #[ignore] // expensive
    fn calling_wham_2d() {;
        Command::new("./target/debug/wham")
            .args(&["--bins", "100,100", "--max", "3.14,3.14", "--min", "-3.14,-3.14", "-T", "300", "--cyclic"])
            .args(&["-f", "example/2d/metadata.dat"])
            .args(&["-o", "/tmp/wham_test_2d.out"])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata("/tmp/wham_test_2d.out").is_ok());
        let output = Command::new("diff")
            .arg("/tmp/wham_test_2d.out")
            .arg("example/2d/wham.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
    }
}