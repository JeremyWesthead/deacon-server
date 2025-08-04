use assert_cmd::Command;
use deacon::MatchThreshold;
use std::fs;
use std::fs::File;
use std::path::Path;
use std::process::Command as StdCommand;
use tempfile::tempdir;

fn create_test_fasta(path: &Path) {
    let fasta_content = ">seq1\nACGTACGTACGT\n>seq2\nGTACGTACGTAC\n";
    fs::write(path, fasta_content).unwrap();
}

fn create_test_fastq(path: &Path) {
    let fastq_content =
        "@seq1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n@seq2\nGTACGTACGTAC\n+\n~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

fn create_test_paired_fastq(path1: &Path, path2: &Path) {
    let fastq_content1 =
        "@read1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n@read2\nGTACGTACGTAC\n+\n~~~~~~~~~~~~\n";
    let fastq_content2 =
        "@read1\nTGCATGCATGCA\n+\n~~~~~~~~~~~~\n@read2\nCATGCATGCATG\n+\n~~~~~~~~~~~~\n";

    fs::write(path1, fastq_content1).unwrap();
    fs::write(path2, fastq_content2).unwrap();
}

fn build_index(fasta_path: &Path, bin_path: &Path) {
    let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
        .arg("index")
        .arg("build")
        .arg(fasta_path)
        .output()
        .expect("Failed to execute command");

    fs::write(bin_path, output.stdout).expect("Failed to write index file");
    assert!(output.status.success(), "Index build command failed");
}

fn create_test_fasta_sc2(path: &Path) {
    let fasta_content =
        ">mn908947.3_0:60\nATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n";
    fs::write(path, fasta_content).unwrap();
}

fn create_test_fastq_sc2_fwd(path: &Path) {
    let fastq_content = "@mn908947.3_0:60_fwd\n\
        ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

fn create_test_fastq_sc2_rev(path: &Path) {
    let fastq_content = "@mn908947.3_0:60_rev\n\
        AGATCTACAAGAGATCGAAAGTTGGTTGGTTTGTTACCTGGGAAGGTATAAACCTTTAAT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

fn create_test_paired_fastq_sc2_fwd(path1: &Path, path2: &Path) {
    let fastq_content1 = "@mn908947.3_0:60_fwd\n\
        ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    let fastq_content2 = "@mn908947.3_60:120_fwd\n\
        GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path1, fastq_content1).unwrap();
    fs::write(path2, fastq_content2).unwrap();
}

fn create_test_paired_fastq_sc2_rev(path1: &Path, path2: &Path) {
    let fastq_content1 = "@mn908947.3_0:60_rev\n\
        AGATCTACAAGAGATCGAAAGTTGGTTGGTTTGTTACCTGGGAAGGTATAAACCTTTAAT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    let fastq_content2 = "@mn908947.3_60:120_rev\n\
        AGTGCACTAAGCATGCAGCCGAGTGACAGCCACACAGATTTTAAAGTTCGTTTAGAGAAC\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path1, fastq_content1).unwrap();
    fs::write(path2, fastq_content2).unwrap();
}

#[test]
fn test_filter_to_file() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let summary_path = temp_dir.path().join("summary.json");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Run filtering command - sequences too short for k=31, so no matches, will be filtered out
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .arg("--summary")
        .arg(&summary_path)
        .assert()
        .success();

    // Check output and report creation
    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    // With new default behavior: sequences without matches are filtered out (sequences too short for k=31)
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.is_empty(),
        "Output file should be empty - sequences too short for minimizers"
    );
}

#[test]
fn test_filter_to_file_gzip() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.gz");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Check gzipped output file creation
    assert!(output_path.exists(), "Gzipped output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "Gzipped output file is empty"
    );
}

#[test]
fn test_filter_to_file_zstd() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.zst");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Check that zstd output file was created
    assert!(output_path.exists(), "Zstd output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "Zstd output file is empty"
    );
}

#[test]
fn test_filter_to_file_xz() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.xz");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Check that xz output file was created
    assert!(output_path.exists(), "XZ output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "XZ output file is empty"
    );
}

#[test]
fn test_filter_deplete_flag() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_depleted.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with deplete flag wasn't created"
    );
}

#[test]
fn test_filter_rename() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_renamed.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--rename")
        .arg("-m")
        .arg("0")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with rename flag wasn't created"
    );

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.contains("@1\n") || output_content.contains("@2\n"),
        "Output does not contain renamed sequences"
    );
}

#[test]
fn test_filter_min_matches() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_min_matches.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("2")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with min_matches parameter wasn't created"
    );
}

#[test]
fn test_filter_prefix_length() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_prefix.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--prefix-length")
        .arg("6")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with prefix_length parameter wasn't created"
    );
}

#[test]
fn test_filter_paired() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Run filtering command with paired-end reads (using -m 0 so short sequences pass through)
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("-m")
        .arg("0")
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Check output creation
    assert!(output_path.exists(), "Output file wasn't created");

    // Validate output content (should be interleaved)
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(!output_content.is_empty(), "Output file is empty");
}

#[test]
fn test_filter_paired_with_deplete() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_depleted.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with deplete flag wasn't created"
    );
}

#[test]
fn test_filter_paired_with_rename() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_renamed.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--rename")
        .arg("-m")
        .arg("0")
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with rename flag wasn't created"
    );

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.contains("@1\n") && output_content.contains("@2\n"),
        "Output does not contain renamed sequences"
    );
}

#[test]
fn test_filter_paired_with_min_matches() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_min_matches.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("2")
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with min_matches parameter wasn't created"
    );
}

#[test]
fn test_interleaved_paired_reads() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let interleaved_fastq_path = temp_dir.path().join("interleaved_reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    // Create test files
    create_test_fasta(&fasta_path);

    let interleaved_content =
        "@read1/1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n@read1/2\nTGCATGCATGCA\n+\n~~~~~~~~~~~~\n"
            .to_owned()
            + "@read2/1\nGTACGTACGTAC\n+\n~~~~~~~~~~~~\n@read2/2\nCATGCATGCATG\n+\n~~~~~~~~~~~~\n";
    fs::write(&interleaved_fastq_path, interleaved_content).unwrap();

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Test piping interleaved file to stdin for processing
    let mut cmd = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"));
    let output = cmd
        .arg("filter")
        .arg("-m")
        .arg("0")
        .arg(&bin_path)
        .arg("-") // stdin for input
        .arg("-") // stdin for input2 (signals interleaved mode)
        .arg("--output")
        .arg(&output_path)
        .stdin(File::open(&interleaved_fastq_path).unwrap())
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command failed");
    assert!(output_path.exists(), "Output file wasn't created");

    // Validate output content (should contain processed reads)
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(!output_content.is_empty(), "Output file is empty");
}

#[test]
fn test_filter_filtration_fwd() {
    // Tests filtering with forward reads from SC2
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let summary_path = temp_dir.path().join("summary.json");

    create_test_fasta_sc2(&fasta_path);
    create_test_fastq_sc2_fwd(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .arg("--summary")
        .arg(&summary_path)
        .arg("--matches")
        .arg("1")
        .assert()
        .success();

    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[test]
fn test_filter_filtration_rev() {
    // Tests filtering with reverse read from SC2
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let summary_path = temp_dir.path().join("summary.json");

    create_test_fasta_sc2(&fasta_path);
    create_test_fastq_sc2_rev(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .arg("--summary")
        .arg(&summary_path)
        .assert()
        .success();

    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[test]
fn test_filter_paired_filtration_fwd() {
    // Tests that both reads are filtered when a forward read matches the SC2 ref
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta_sc2(&fasta_path);
    create_test_paired_fastq_sc2_fwd(&fastq_path1, &fastq_path2);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(output_path.exists(), "Output file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[test]
fn test_filter_paired_filtration_rev() {
    // Tests that both reads are filtered when a reverse read matches the SC2 ref
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta_sc2(&fasta_path);
    create_test_paired_fastq_sc2_rev(&fastq_path1, &fastq_path2);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(output_path.exists(), "Output file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[cfg(test)]
mod output2_tests {
    use assert_cmd::Command;
    use std::fs;
    use std::path::Path;
    use std::process::Command as StdCommand;
    use tempfile::tempdir;

    fn create_test_fasta(path: &Path) {
        let fasta_content = ">seq1\nACGTACGTACGT\n>seq2\nGTACGTACGTAC\n";
        fs::write(path, fasta_content).unwrap();
    }

    fn create_test_paired_fastq(path1: &Path, path2: &Path) {
        let fastq_content1 =
            "@read1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n@read2\nGTACGTACGTAC\n+\n~~~~~~~~~~~~\n";
        let fastq_content2 =
            "@read1\nTGCATGCATGCA\n+\n~~~~~~~~~~~~\n@read2\nCATGCATGCATG\n+\n~~~~~~~~~~~~\n";

        fs::write(path1, fastq_content1).unwrap();
        fs::write(path2, fastq_content2).unwrap();
    }

    fn build_index(fasta_path: &Path, bin_path: &Path) {
        let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
            .arg("index")
            .arg("build")
            .arg(fasta_path)
            .output()
            .expect("Failed to execute command");

        fs::write(bin_path, output.stdout).expect("Failed to write index file");
        assert!(output.status.success(), "Index build command failed");
    }

    #[test]
    fn test_filter_paired_with_output2() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("ref.fasta");
        let fastq_path1 = temp_dir.path().join("reads_1.fastq");
        let fastq_path2 = temp_dir.path().join("reads_2.fastq");
        let bin_path = temp_dir.path().join("ref.bin");
        let output_path1 = temp_dir.path().join("filtered_1.fastq");
        let output_path2 = temp_dir.path().join("filtered_2.fastq");
        let summary_path = temp_dir.path().join("summary.json");

        create_test_fasta(&fasta_path);
        create_test_paired_fastq(&fastq_path1, &fastq_path2);

        build_index(&fasta_path, &bin_path);
        assert!(bin_path.exists(), "Index file wasn't created");

        // Run filtering command with separate output files (using -m 0 so short sequences pass through)
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("filter")
            .arg("-m")
            .arg("0")
            .arg(&bin_path)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path1)
            .arg("--output2")
            .arg(&output_path2)
            .arg("--summary")
            .arg(&summary_path)
            .assert()
            .success();

        // Check both output files were created
        assert!(output_path1.exists(), "First output file wasn't created");
        assert!(output_path2.exists(), "Second output file wasn't created");
        assert!(summary_path.exists(), "Summary file wasn't created");

        // Validate output content
        let output1_content = fs::read_to_string(&output_path1).unwrap();
        let output2_content = fs::read_to_string(&output_path2).unwrap();

        assert!(!output1_content.is_empty(), "First output file is empty");
        assert!(!output2_content.is_empty(), "Second output file is empty");

        // Check that the summary includes output2 path
        let summary_content = fs::read_to_string(&summary_path).unwrap();
        assert!(
            summary_content.contains("output2"),
            "Summary doesn't mention output2"
        );
    }

    #[test]
    fn test_filter_paired_with_output2_gzip() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("ref.fasta");
        let fastq_path1 = temp_dir.path().join("reads_1.fastq");
        let fastq_path2 = temp_dir.path().join("reads_2.fastq");
        let bin_path = temp_dir.path().join("ref.bin");
        let output_path1 = temp_dir.path().join("filtered_1.fastq.gz");
        let output_path2 = temp_dir.path().join("filtered_2.fastq.gz");

        create_test_fasta(&fasta_path);
        create_test_paired_fastq(&fastq_path1, &fastq_path2);
        build_index(&fasta_path, &bin_path);

        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("filter")
            .arg(&bin_path)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path1)
            .arg("--output2")
            .arg(&output_path2)
            .assert()
            .success();

        // Check both gzipped output files were created
        assert!(
            output_path1.exists(),
            "First gzipped output file wasn't created"
        );
        assert!(
            output_path2.exists(),
            "Second gzipped output file wasn't created"
        );

        assert!(
            fs::metadata(&output_path1).unwrap().len() > 0,
            "First gzipped output file is empty"
        );
        assert!(
            fs::metadata(&output_path2).unwrap().len() > 0,
            "Second gzipped output file is empty"
        );
    }

    #[test]
    fn test_filter_single_input_with_output2_warning() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("ref.fasta");
        let fastq_path = temp_dir.path().join("reads.fastq");
        let bin_path = temp_dir.path().join("ref.bin");
        let output_path = temp_dir.path().join("filtered.fastq");
        let output_path2 = temp_dir.path().join("filtered_2.fastq");

        create_test_fasta(&fasta_path);
        let fastq_content = "@seq1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n";
        fs::write(&fastq_path, fastq_content).unwrap();
        build_index(&fasta_path, &bin_path);

        // Run filtering command with output2 but no second input (should warn)
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("filter")
            .arg(&bin_path)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .arg("-O")
            .arg(&output_path2)
            .assert()
            .success()
            .stderr(predicates::str::contains("Warning"));

        // Check only the first output file was created
        assert!(output_path.exists(), "First output file wasn't created");
        assert!(
            !output_path2.exists(),
            "Second output file shouldn't be created for single input"
        );
    }
}

#[test]
fn test_shared_minimizer_counted_once() {
    // Catch bug where the same minimizer in different paired mates is counted twice
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fasta_path1 = temp_dir.path().join("reads_1.fasta");
    let fasta_path2 = temp_dir.path().join("reads_2.fasta");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fasta");
    let summary_path = temp_dir.path().join("summary.json");

    let ref_content = ">reference\nACGTACGTACGTACGTTGCATGCATGCATGCATAAGGTTAAGGTTAAGGTTAAGGTTCCCGGGCCCGGGCCCGGGCCCGGGATATATATATATATATATGCGCGCGCGCGCGCGCGC\n";
    fs::write(&fasta_path, ref_content).unwrap();

    // Create 120bp ref
    let ref_content = ">reference\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    fs::write(&fasta_path, ref_content).unwrap();

    // Create paired reads (80bp each) where both contain the same 60bp region from the reference
    // Shared region: first 60bp of reference (ACGT repeated 15 times)
    let fasta_content1 = ">read1/1\n\
        AAAAAAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAAAAAAA\n";

    let fasta_content2 = ">read1/2\n\
        TTTTTTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTT\n";

    fs::write(&fasta_path1, fasta_content1).unwrap();
    fs::write(&fasta_path2, fasta_content2).unwrap();

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // If shared minimizers are counted once (correct): total hits = 1, pair kept (1 < 2)
    // If shared minimizers are counted twice (bug): total hits = 2+, pair filtered (2+ >= 2)
    // Using --deplete to restore original behavior for this bug test
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--deplete")
        .arg(&bin_path)
        .arg(&fasta_path1)
        .arg(&fasta_path2)
        .arg("--output")
        .arg(&output_path)
        .arg("--summary")
        .arg(&summary_path)
        .arg("--matches")
        .arg("2") // Critical parameter: any pair with 2+ hits gets filtered
        .assert()
        .success();

    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    let summary_content = fs::read_to_string(&summary_path).unwrap();
    let summary: serde_json::Value = serde_json::from_str(&summary_content).unwrap();

    // The reads should be kept because shared minimizers should only count once
    assert!(
        !output_content.is_empty(),
        "Read pair should be kept in output because shared minimizers should only count once. \
         Current implementation incorrectly counts them multiple times and filters the pair."
    );

    // Additional verification using the JSON summary
    let seqs_out = summary["seqs_out"].as_u64().unwrap();
    assert_eq!(
        seqs_out, 2,
        "Expected 2 sequences in output (both reads of the pair should be kept) \
         but got {seqs_out}. This indicates shared minimizers were double-counted."
    );
}

#[test]
fn test_match_threshold_parsing() {
    // Test MatchThreshold parsing functionality
    assert_eq!(
        "2".parse::<MatchThreshold>().unwrap(),
        MatchThreshold::Absolute(2)
    );
    assert_eq!(
        "0.5".parse::<MatchThreshold>().unwrap(),
        MatchThreshold::Relative(0.5)
    );
    assert_eq!(
        "0.0".parse::<MatchThreshold>().unwrap(),
        MatchThreshold::Relative(0.0)
    );
    assert_eq!(
        "1.0".parse::<MatchThreshold>().unwrap(),
        MatchThreshold::Relative(1.0)
    );

    // Test invalid inputs
    assert!("1.5".parse::<MatchThreshold>().is_err()); // > 1.0
    assert!("-0.1".parse::<MatchThreshold>().is_err()); // < 0.0
    assert!("abc".parse::<MatchThreshold>().is_err()); // not a number
}

#[test]
fn test_match_threshold_logic() {
    let abs_threshold = MatchThreshold::Absolute(3);
    let prop_threshold = MatchThreshold::Relative(0.5);

    // Test absolute threshold
    let required_abs = match abs_threshold {
        MatchThreshold::Absolute(n) => n,
        _ => panic!(),
    };
    assert!(3 >= required_abs);
    assert!(5 >= required_abs);
    assert!(2 < required_abs);

    // Test relative threshold (50% of 10 minimizers = 5 hits required)
    let required_rel = match prop_threshold {
        MatchThreshold::Relative(f) => ((f * 10.0).ceil() as usize).max(1),
        _ => panic!(),
    };
    assert!(5 >= required_rel);
    assert!(6 >= required_rel);
    assert!(4 < required_rel);

    // Test edge case: minimum 1 hit required even for small proportions
    let small_prop = MatchThreshold::Relative(0.1);
    let required_small = match small_prop {
        MatchThreshold::Relative(f) => ((f * 5.0).ceil() as usize).max(1),
        _ => panic!(),
    };
    assert!(1 >= required_small); // 0.1 * 5 = 0.5, ceil to 1
    assert!(0 < required_small);

    // Test edge case: zero minimizers case would be handled at runtime
}

#[test]
fn test_filter_proportional_threshold() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_proportional.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("0.5") // 50% proportional threshold
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with proportional threshold wasn't created"
    );
}

#[test]
fn test_filter_proportional_paired() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_proportional_paired.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("0.3") // 30% proportional threshold
        .arg(&bin_path)
        .arg(&fastq_path1)
        .arg(&fastq_path2)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with proportional threshold for paired reads wasn't created"
    );
}

#[test]
fn test_filter_edge_case_proportional_values() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_edge.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    // Test with 0.0 (should pass everything)
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("0.0")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Test with 1.0 (very strict)
    let output_path_strict = temp_dir.path().join("filtered_strict.fastq");
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("1.0")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path_strict)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output with 0.0 threshold wasn't created"
    );
    assert!(
        output_path_strict.exists(),
        "Output with 1.0 threshold wasn't created"
    );
}
