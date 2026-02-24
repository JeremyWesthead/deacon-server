use assert_cmd::cargo;
use predicates::str;

#[test]
fn test_version() {
    let mut cmd = cargo::cargo_bin_cmd!("deacon");
    cmd.arg("--version")
        .assert()
        .success()
        .stdout(str::contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn test_no_args() {
    let mut cmd = cargo::cargo_bin_cmd!("deacon");
    cmd.assert().failure().stderr(str::contains("Usage"));
}
