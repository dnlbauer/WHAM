use std::process::Command;
use std::env;

pub fn get_command() -> Command {
    let exe = env::current_exe().unwrap();
    let path = exe.to_str().unwrap();
    let cmd = if path.contains("release") {
        "./target/release/wham"
    } else {
        "./target/debug/wham"
    };
    Command::new(cmd)
}