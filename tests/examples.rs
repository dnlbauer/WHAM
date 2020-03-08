mod command;

#[cfg(test)]
mod integration {

    use std::process::Command;
    use std::fs;
    use super::command::get_command;

    #[test]
    fn wham_1d() {
        get_command()
            .args(&["--bins", "100", "--max", "4.5", "--min", "0.7", "-T", "310"])
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
    fn wham_1d_cyclic() {
        get_command()
            .args(&["--bins", "100", "--max", "pi", "--min", "-pi", "-T", "300", "--cyclic"])
            .args(&["-f", "example/1d_cyclic/metadata.dat"])
            .args(&["-o", "/tmp/wham_test_1d_cyclic.out"])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata("/tmp/wham_test_1d_cyclic.out").is_ok());
        let output = Command::new("diff")
            .arg("/tmp/wham_test_1d_cyclic.out")
            .arg("example/1d_cyclic/wham.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
    }

    #[test]
    #[ignore] // expensive
    fn wham_2d_cyclic() {
        get_command()
            .args(&["--bins", "100,100", "--max", "pi,pi", "--min", "-pi,-pi", "-T", "300", "--cyclic"])
            .args(&["-f", "example/2d_cyclic/metadata.dat"])
            .args(&["-o", "/tmp/wham_test_2d_cyclic.out"])
            .output()
            .expect("failed to execute process");

        assert!(fs::metadata("/tmp/wham_test_2d_cyclic.out").is_ok());
        let output = Command::new("diff")
            .arg("/tmp/wham_test_2d_cyclic.out")
            .arg("example/2d_cyclic/wham.out")
            .output()
            .expect("failed to run diff");
        let output_len = String::from_utf8_lossy(&output.stdout).len();
        assert_eq!(output_len, 0);
    }
}