#[cfg(feature = "server")]
use crate::MatchThreshold;
use crate::index::load_minimizer_hashes;
use crate::minimizers::fill_minimizer_hashes;
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use liblzma::write::XzEncoder;
use needletail::parse_fastx_file;
use needletail::parse_fastx_stdin;
use needletail::parser::Format;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use zstd::stream::write::Encoder as ZstdEncoder;

#[cfg(feature = "server")]
use crate::server_common::{FilterRequest, FilterResponse};
#[cfg(feature = "server")]
use reqwest::Client;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer

struct RecordData {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
    format: Format,
}
trait FastxWriter: Write {
    fn flush_all(&mut self) -> io::Result<()>;
}

trait CompressionEncoder: Write {
    fn finish(self: Box<Self>) -> io::Result<()>;
}

#[derive(Debug, Clone, Copy)]
enum CompressionFormat {
    None,
    Gzip,
    Zstd,
    Xz,
}

impl CompressionFormat {
    fn from_extension(path: &str) -> Self {
        if path.ends_with(".gz") {
            Self::Gzip
        } else if path.ends_with(".zst") {
            Self::Zstd
        } else if path.ends_with(".xz") {
            Self::Xz
        } else {
            Self::None
        }
    }

    fn validate_compression_level(&self, level: u8) -> Result<()> {
        match self {
            Self::None => Ok(()),
            Self::Gzip => {
                if level < 1 || level > 9 {
                    Err(anyhow::anyhow!(
                        "Invalid gzip compression level {}. Must be between 1 and 9.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
            Self::Zstd => {
                if level < 1 || level > 22 {
                    Err(anyhow::anyhow!(
                        "Invalid zstd compression level {}. Must be between 1 and 22.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
            Self::Xz => {
                if level > 9 {
                    Err(anyhow::anyhow!(
                        "Invalid xz compression level {}. Must be between 0 and 9.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }
}

impl<W: Write> CompressionEncoder for GzEncoder<W> {
    fn finish(mut self: Box<Self>) -> io::Result<()> {
        self.try_finish()
    }
}

impl<W: Write> CompressionEncoder for ZstdEncoder<'static, W> {
    fn finish(self: Box<Self>) -> io::Result<()> {
        (*self).finish().map(|_| ())
    }
}

impl<W: Write> CompressionEncoder for XzEncoder<W> {
    fn finish(self: Box<Self>) -> io::Result<()> {
        (*self).finish().map(|_| ())
    }
}

struct CompressedWriter {
    encoder: Option<Box<dyn CompressionEncoder>>,
}

impl CompressedWriter {
    fn new(encoder: Box<dyn CompressionEncoder>) -> Self {
        Self {
            encoder: Some(encoder),
        }
    }

    fn uncompressed<W: Write>(writer: W) -> StandardWriter<W> {
        StandardWriter(writer)
    }
}

impl Write for CompressedWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Some(encoder) = &mut self.encoder {
            encoder.write(buf)
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been closed",
            ))
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(encoder) = &mut self.encoder {
            encoder.flush()
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been closed",
            ))
        }
    }
}

impl FastxWriter for CompressedWriter {
    fn flush_all(&mut self) -> io::Result<()> {
        if let Some(encoder) = self.encoder.take() {
            encoder.finish()?;
        }
        Ok(())
    }
}

struct StandardWriter<W: Write>(W);

impl<W: Write> Write for StandardWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

impl<W: Write> FastxWriter for StandardWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        self.flush()
    }
}

// Return a file writer appropriate for the output path extension
fn get_writer(output_path: &str, compression_level: u8) -> Result<Box<dyn FastxWriter>> {
    if output_path == "-" {
        // Write to stdout
        let stdout = io::stdout();
        let writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, stdout);
        Ok(Box::new(CompressedWriter::uncompressed(writer)))
    } else {
        // Write to file with extension-appropriate encoder
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(output_path)
            .context(format!("Failed to create output file: {}", output_path))?;

        let buffered_file = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);
        let format = CompressionFormat::from_extension(output_path);

        // Validate compression level for the format
        format.validate_compression_level(compression_level)?;

        match format {
            CompressionFormat::None => Ok(Box::new(CompressedWriter::uncompressed(buffered_file))),
            CompressionFormat::Gzip => {
                let encoder =
                    GzEncoder::new(buffered_file, Compression::new(compression_level as u32));
                Ok(Box::new(CompressedWriter::new(Box::new(encoder))))
            }
            CompressionFormat::Zstd => {
                let encoder = ZstdEncoder::new(buffered_file, compression_level as i32)
                    .context("Failed to create zstd encoder")?;
                Ok(Box::new(CompressedWriter::new(Box::new(encoder))))
            }
            CompressionFormat::Xz => {
                let encoder = XzEncoder::new(buffered_file, compression_level as u32);
                Ok(Box::new(CompressedWriter::new(Box::new(encoder))))
            }
        }
    }
}

// JSON summary structure
#[derive(Serialize, Deserialize)]
pub struct FilterSummary {
    version: String,
    index: String,
    input: String,
    input2: Option<String>,
    output: String,
    output2: Option<String>,
    k: usize,
    w: usize,
    match_threshold: String,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    seqs_in: u64,
    seqs_out: u64,
    seqs_out_proportion: f64,
    seqs_removed: u64,
    seqs_removed_proportion: f64,
    bp_in: u64,
    bp_out: u64,
    bp_out_proportion: f64,
    bp_removed: u64,
    bp_removed_proportion: f64,
    time: f64,
    seqs_per_second: u64,
    bp_per_second: u64,
}

/// Determine if a given set of minimizers matches the index
pub fn input_matches_index(
    index_minimizers: &FxHashSet<u64>,
    input_minimizers: &Vec<u64>,
    match_threshold: &MatchThreshold,
) -> bool {
    // Count distinct minimizer hits
    let mut seen_hits = FxHashSet::default();
    let mut hit_count = 0;
    for &hash in input_minimizers {
        if index_minimizers.contains(&hash) && seen_hits.insert(hash) {
            hit_count += 1;
        }
    }

    // Convert threshold to absolute count
    let required_hits = match match_threshold {
        MatchThreshold::Absolute(n) => *n,
        MatchThreshold::Relative(f) => {
            if input_minimizers.is_empty() {
                0
            } else {
                ((*f * input_minimizers.len() as f64).ceil() as usize).max(1)
            }
        }
    };

    // The input minimizers match the index if it has enough hits
    hit_count >= required_hits
}

pub fn inputs_match_index(
    index_minimizers: &FxHashSet<u64>,
    input_minimizers: &Vec<Vec<u64>>,
    match_threshold: &MatchThreshold,
) -> Vec<bool> {
    input_minimizers
        .par_iter()
        .map(|minimizers| {
            input_matches_index(index_minimizers, minimizers, match_threshold)
        })
        .collect()
}

/// Send minimizers to server for checking against index.
/// Equivalent functionality to `input_matches_index, but remote
#[cfg(feature = "server")]
async fn send_minimizers_to_server(
    input_minimizers: Vec<u64>,
    server_address: &str,
    matches_threshold: &MatchThreshold,
) -> Result<bool> {
    // Create a client to send the minimizers to the server
    let client = Client::new();

    // Send the minimizers as a POST request
    let response = client
        .post(server_address.to_owned() + "/is_index_match")
        .json(&FilterRequest {
            input: vec![input_minimizers],
            match_threshold: *matches_threshold,
        })
        .send()
        .await?;

    // Check if the response indicates a match
    if response.status().is_success() {
        Ok(response.json::<FilterResponse>().await?.index_match[0])
    } else {
        Err(anyhow::anyhow!(
            "Server returned an error: {}",
            response.status()
        ))
    }
}

/// Send minimizers to server for checking against index.
/// Equivalent functionality to `input_matches_index, but remote
#[cfg(feature = "server")]
async fn send_all_minimizers_to_server(
    input_minimizers: Vec<Vec<u64>>,
    server_address: &str,
    matches_threshold: &MatchThreshold,
) -> Result<Vec<bool>> {
    // Create a client to send the minimizers to the server
    let client = Client::new();

    // Send the minimizers as a POST request
    let response = client
        .post(server_address.to_owned() + "/is_index_match")
        .json(&FilterRequest {
            input: input_minimizers,
            match_threshold: *matches_threshold,
        })
        .send()
        .await?;

    // Check if the response indicates a match
    if response.status().is_success() {
        Ok(response.json::<FilterResponse>().await?.index_match)
    } else {
        Err(anyhow::anyhow!(
            "Server returned an error: {}",
            response.status()
        ))
    }
}

pub async fn check_input_matches_index(
    index_minimizers: &Option<FxHashSet<u64>>,
    input_minimizers: &Vec<u64>,
    match_threshold: &MatchThreshold,
    server_address: &Option<String>,
) -> bool {
    // If index minimizers are provided, check if input matches locally
    if let Some(index_minimizers) = index_minimizers {
        input_matches_index(index_minimizers, input_minimizers, match_threshold)
    } else {
        // Else, send the input minimizers to the server for checking
        #[cfg(feature = "server")]
        {
            if server_address.is_none() {
                panic!("Server address is required when using the server feature.");
            }
            let server_address = server_address.as_ref().map(String::as_str).unwrap();
            return send_minimizers_to_server(
                input_minimizers.to_vec(),
                server_address,
                match_threshold,
            )
            .await
            .unwrap_or_else(|e| {
                panic!("Error checking input against index: {}", e);
            });
        }
        #[cfg(not(feature = "server"))]
        {
            panic!("Server feature is not enabled. Cannot check input against index.");
        }
    }
}

//TODO: Update to use just this function
// Sever appears to be blocking when the index minimizers are accessed :(
// So group together the batches of input minimizers and send them all at once
// Multithreaded server access means we get a resonable filter speed
pub async fn check_inputs_match_index(
    index_minimizers: &Option<FxHashSet<u64>>,
    input_minimizers: &Vec<Vec<u64>>,
    match_threshold: &MatchThreshold,
    server_address: &Option<String>,
) -> Vec<bool> {
    // If index minimizers are provided, check if input matches locally
    if let Some(index_minimizers) = index_minimizers {
        inputs_match_index(index_minimizers, input_minimizers, match_threshold)
    } else {
        // Else, send the input minimizers to the server for checking
        #[cfg(feature = "server")]
        {
            if server_address.is_none() {
                panic!("Server address is required when using the server feature.");
            }
            let server_address = server_address.as_ref().map(String::as_str).unwrap();
            return send_all_minimizers_to_server(
                input_minimizers.to_vec(),
                server_address,
                match_threshold,
            )
            .await
            .unwrap_or_else(|e| {
                panic!("Error checking input against index: {}", e);
            });
        }
        #[cfg(not(feature = "server"))]
        {
            panic!("Server feature is not enabled. Cannot check input against index.");
        }
    }
}

pub async fn run<P: AsRef<Path>>(
    minimizers_path: Option<P>,
    input_path: &str,
    input2_path: Option<&str>,
    output_path: &str,
    output2_path: Option<&str>,
    match_threshold: &MatchThreshold,
    prefix_length: usize,
    summary_path: Option<&PathBuf>,
    deplete: bool,
    rename: bool,
    threads: usize,
    compression_level: u8,
    server_address: Option<String>,
) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    let tool_version = format!("deacon {}", version);

    // Configure thread pool if nonzero
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to initialize thread pool")?;
    }

    let mode = if deplete { "deplete" } else { "search" };

    let mut input_type = String::new();
    let mut options = Vec::<String>::new();
    let paired_stdin = input_path == "-" && input2_path.is_some() && input2_path.unwrap() == "-";
    if paired_stdin {
        input_type.push_str("interleaved");
    } else if let Some(_) = input2_path {
        input_type.push_str("paired");
    } else {
        input_type.push_str("single");
    }
    options.push(format!("match_threshold={}", match_threshold));
    if prefix_length > 0 {
        options.push(format!("prefix_length={}", prefix_length));
    }
    if rename {
        options.push("rename".to_string());
    }
    if threads > 0 {
        options.push(format!("threads={}", threads));
    }

    eprintln!(
        "Deacon v{}; mode: {}; input: {}; options: {}",
        version,
        mode,
        input_type,
        options.join(", ")
    );

    // Load minimizers hashes and parse header
    let (minimizer_hashes, header) =
        load_minimizer_hashes(&minimizers_path, &server_address).await?;

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let load_time = start_time.elapsed();
    eprintln!(
        "Loaded index (k={}, w={}) in {:.2?}",
        kmer_length, window_size, load_time
    );

    // Create the appropriate writer(s) based on the output path(s)
    let mut writer = get_writer(output_path, compression_level)?;
    let mut writer2 = if let (Some(output2), Some(_)) = (output2_path, input2_path) {
        // Only create second writer if both output2 and input2 are specified
        Some(get_writer(output2, compression_level)?)
    } else {
        None
    };

    // A progress bar would require a denominator, so let's spin
    let spinner = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[".  ", ".. ", "...", " ..", "  .", "   "])
            .template("{msg}{spinner} ")?,
    );
    spinner.set_message("Filtering");

    // Init counters
    let mut total_seqs = 0;
    let mut filtered_seqs = 0;
    let mut total_bp = 0;
    let mut output_bp = 0;
    let mut filtered_bp = 0;
    let mut output_seq_counter = 0;

    // Start timer for filtering rate calculation (excludes index loading time)
    let filtering_start_time = Instant::now();

    if paired_stdin {
        process_interleaved_paired_seqs(
            &minimizer_hashes,
            &mut writer,
            writer2.as_mut(),
            match_threshold,
            prefix_length,
            kmer_length,
            window_size,
            deplete,
            rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            filtering_start_time,
            &server_address,
        )
        .await?;
    } else if let Some(input2_path) = input2_path {
        process_paired_seqs(
            &minimizer_hashes,
            input_path,
            input2_path,
            &mut writer,
            writer2.as_mut(),
            match_threshold,
            prefix_length,
            kmer_length,
            window_size,
            deplete,
            rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            filtering_start_time,
            &server_address,
        )
        .await?;
    } else {
        process_single_seqs(
            &minimizer_hashes,
            input_path,
            &mut writer,
            match_threshold,
            prefix_length,
            kmer_length,
            window_size,
            deplete,
            rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            filtering_start_time,
            &server_address,
        )
        .await?;
    }

    writer.flush_all()?;
    if let Some(ref mut w2) = writer2 {
        w2.flush_all()?;
    }

    let total_time = start_time.elapsed();
    let seqs_per_sec = total_seqs as f64 / total_time.as_secs_f64();
    let bp_per_sec = total_bp as f64 / total_time.as_secs_f64();
    let mbp_per_sec = bp_per_sec / 1_000_000.0;

    // Calculate filtered proportion directly
    let filtered_proportion = if total_seqs > 0 {
        filtered_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    // Calculate filtered base pair proportion
    let filtered_bp_proportion = if total_bp > 0 {
        filtered_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    // Calculate output proportions
    let output_seqs = total_seqs - filtered_seqs;
    let output_seq_proportion = if total_seqs > 0 {
        output_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    let output_bp_proportion = if total_bp > 0 {
        output_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    // Finish and clear spinner, print final message
    spinner.finish_and_clear();
    eprintln!(
        "Retained {}/{} sequences ({:.3}%), {}/{} bp ({:.3}%)",
        output_seqs,
        total_seqs,
        output_seq_proportion * 100.0,
        output_bp,
        total_bp,
        output_bp_proportion * 100.0
    );

    // Print completion message with speed
    eprintln!(
        "Completed in {:.2?}. Speed: {:.0} seqs/s ({:.1} Mbp/s)",
        total_time, seqs_per_sec, mbp_per_sec
    );

    // Build and write a JSON summary if path provided
    if let Some(summary_file) = summary_path {
        // Get number of sequences passing filter
        let seqs_out = total_seqs - filtered_seqs;

        let index = match minimizers_path {
            Some(path) => path.as_ref().to_string_lossy().to_string(),
            None => server_address.unwrap(),
        };

        let summary = FilterSummary {
            version: tool_version,
            index,
            input: input_path.to_string(),
            input2: input2_path.map(|s| s.to_string()),
            output: output_path.to_string(),
            output2: output2_path.map(|s| s.to_string()),
            k: kmer_length,
            w: window_size,
            match_threshold: match_threshold.to_string(),
            prefix_length,
            deplete,
            rename,
            seqs_in: total_seqs as u64,
            seqs_out: seqs_out as u64,
            seqs_out_proportion: output_seq_proportion,
            seqs_removed: filtered_seqs as u64,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp as u64,
            bp_out: output_bp as u64,
            bp_out_proportion: output_bp_proportion,
            bp_removed: filtered_bp as u64,
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
        };

        // Write summary file
        let file = File::create(summary_file)
            .context(format!("Failed to create summary: {:?}", summary_file))?;
        let writer = BufWriter::new(file);

        // Serialise and write the summary JSON
        serde_json::to_writer_pretty(writer, &summary).context("Failed to write summary")?;

        eprintln!("Summary saved to {:?}", summary_file);
    }

    Ok(())
}

async fn record_match(
    record_hashes: &Vec<u64>,
    minimizer_hashes: &Option<FxHashSet<u64>>,
    match_threshold: &MatchThreshold,
    server_address: &Option<String>,
    deplete: bool,
) -> bool {
    // check_input_matches_index returns true if the matches >= match_threshold
    // If deplete is true, we want to suppress output of sequences that match
    // The NOT is required due to the logic compared to requirement
    // matches_index && deplete = true --> false
    // !matches_index && deplete = false --> true
    // matches_index && !deplete = false --> true
    // !matches_index && !deplete = true --> false
    let should_output = !(check_input_matches_index(
        minimizer_hashes,
        record_hashes,
        match_threshold,
        server_address,
    )
    .await
        && deplete);

    should_output
}

fn get_hashes_from_record(
    record_data: &RecordData,
    kmer_length: usize,
    prefix_length: usize,
    window_size: usize,
) -> (Vec<u64>, usize) {
    let seq_len = record_data.seq.len();
    let mut minimizer_buffer = Vec::with_capacity(64);

    if seq_len >= kmer_length {
        // Apply prefix length limit if specified
        let effective_seq = if prefix_length > 0 && seq_len > prefix_length {
            &record_data.seq[..prefix_length]
        } else {
            &record_data.seq
        };

        // Get minimizer hash values using parameters from header
        fill_minimizer_hashes(
            effective_seq,
            kmer_length,
            window_size,
            &mut minimizer_buffer,
        );
    }

    (minimizer_buffer, seq_len)
}

fn separate_tuple_vec<T, T2>(input_vec: Vec<(T, T2)>) -> (Vec<T>, Vec<T2>) {
    let mut first = Vec::with_capacity(input_vec.len());
    let mut second = Vec::with_capacity(input_vec.len());

    for (a, b) in input_vec {
        first.push(a);
        second.push(b);
    }

    (first, second)
}

async fn process_single_seqs(
    minimizer_hashes: &Option<FxHashSet<u64>>,
    input_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    match_threshold: &MatchThreshold,
    prefix_length: usize,
    kmer_length: usize,
    window_size: usize,
    deplete: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    filtering_start_time: Instant,
    server_address: &Option<String>,
) -> Result<()> {
    // Create a reader based on the input source
    let mut reader = if input_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input_path)?
    };

    // Process in batches
    let batch_size = 10000;
    let mut output_record_buffer = Vec::with_capacity(1024);

    // Process batches
    loop {
        // Collect a batch of records with owned data
        let mut batch: Vec<RecordData> = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        // Fill the batch (sequential read from either stdin or file)
        for _ in 0..batch_size {
            if let Some(record_result) = reader.next() {
                match record_result {
                    Ok(record) => {
                        let record_data = RecordData {
                            id: record.id().to_vec(),
                            seq: record.seq().to_vec(),
                            qual: record.qual().map(|q| q.to_vec()),
                            format: record.format(),
                        };
                        batch.push(record_data);
                    }
                    Err(e) => return Err(e.into()),
                }
            } else {
                reached_end = true;
                break;
            }
        }

        if batch.is_empty() {
            break;
        }

        // Process batch in parallel
        // let batch_results: Vec<_> = batch
        //     .par_iter()
        //     .map(|record_data| {
        //         record_match(
        //             record_data,
        //             kmer_length,
        //             prefix_length,
        //             window_size,
        //             minimizer_hashes,
        //             match_threshold,
        //             server_address,
        //             deplete,
        //         )
        //     })
        //     .collect();

        let batch_result = batch
            .par_iter()
            .map(|record_data| {
                get_hashes_from_record(
                    record_data,
                    kmer_length,
                    prefix_length,
                    window_size,
                )
            })
            .collect::<Vec<_>>();
        
        let (batch_hashes, seq_lens): (Vec<_>, Vec<_>) = separate_tuple_vec(batch_result);
        
        
        let batch_matches = check_inputs_match_index(
            minimizer_hashes,
            &batch_hashes,
            match_threshold,
            server_address,
        ).await;
        
        

        // Process results sequentially to maintain order
        for (i, (seq_len, seq_match)) in seq_lens.iter().zip(batch_matches).enumerate() {
            let record_data = &batch[i];
            let should_output = !(seq_match && deplete);
            *total_seqs += 1;
            *total_bp += *seq_len as u64;

            if should_output {
                // Track output base pairs
                *output_bp += *seq_len as u64;

                // Increment output sequence counter
                *output_seq_counter += 1;

                // Format as FASTX and write
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    &record_data.id,
                    &record_data.seq,
                    record_data.qual.as_deref(),
                    record_data.format,
                    &mut output_record_buffer,
                    rename,
                    *output_seq_counter,
                );
                writer.write_all(&output_record_buffer)?;
            } else {
                *filtered_seqs += 1;
                *filtered_bp += *seq_len as u64;
            }
        }

        // Update spinner and flush periodically
        let elapsed = filtering_start_time.elapsed();
        let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
        let mbp_per_sec = bp_per_sec / 1_000_000.0;

        // Calculate output proportion
        let output_seqs = *total_seqs - *filtered_seqs;
        let output_proportion = if *total_seqs > 0 {
            output_seqs as f64 / *total_seqs as f64
        } else {
            0.0
        };

        // Calculate output base pair proportion
        let output_bp_proportion = if *total_bp > 0 {
            *output_bp as f64 / *total_bp as f64
        } else {
            0.0
        };

        // Update spinner message
        spinner.set_message(format!(
            "Retained {}/{} sequences ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            seqs_per_sec,
            mbp_per_sec
        ));

        // Flush writer periodically
        writer.flush()?;

        // Check if we've reached the end of the file/stdin
        if reached_end {
            break;
        }
    }

    Ok(())
}

async fn process_paired_seqs(
    minimizer_hashes: &Option<FxHashSet<u64>>,
    input1_path: &str,
    input2_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    mut writer2: Option<&mut Box<dyn FastxWriter>>,
    match_threshold: &MatchThreshold,
    prefix_length: usize,
    kmer_length: usize,
    window_size: usize,
    deplete: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    filtering_start_time: Instant,
    server_address: &Option<String>,
) -> Result<()> {
    // Open both input files
    let mut reader1 = if input1_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input1_path)?
    };

    let mut reader2 = parse_fastx_file(input2_path)?;

    // Process in batches
    let batch_size = 10000;
    let mut output_record_buffer = Vec::with_capacity(1024);

    // Process batches
    loop {
        // Collect a batch of read pairs with owned data
        let mut batch1: Vec<RecordData> = Vec::with_capacity(batch_size);
        let mut batch2: Vec<RecordData> = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        // Fill the batch (sequential read from files)
        for _ in 0..batch_size {
            if let (Some(record1_res), Some(record2_res)) = (reader1.next(), reader2.next()) {
                match (record1_res, record2_res) {
                    (Ok(record1), Ok(record2)) => {
                        let record_data1 = RecordData {
                            id: record1.id().to_vec(),
                            seq: record1.seq().to_vec(),
                            qual: record1.qual().map(|q| q.to_vec()),
                            format: record1.format(),
                        };
                        let record_data2 = RecordData {
                            id: record2.id().to_vec(),
                            seq: record2.seq().to_vec(),
                            qual: record2.qual().map(|q| q.to_vec()),
                            format: record2.format(),
                        };
                        batch1.push(record_data1);
                        batch2.push(record_data2);
                    }
                    (Err(e), _) => return Err(e.into()),
                    (_, Err(e)) => return Err(e.into()),
                }
            } else {
                reached_end = true;
                break;
            }
        }

        if batch1.is_empty() {
            break;
        }

        // Process batch in parallel
        let batch_results: Vec<_> = batch1
            .par_iter()
            .zip(batch2.par_iter())
            .map(async |(record_data1, record_data2)| {
                let seq1_len = record_data1.seq.len();
                let seq2_len = record_data2.seq.len();

                // Pre-allocate buffers for reuse
                let mut minimizer_buffer1 = Vec::with_capacity(64);
                let mut minimizer_buffer2 = Vec::with_capacity(64);

                // Check for minimizer hits in read 1
                if seq1_len >= kmer_length {
                    minimizer_buffer1.clear();

                    // Apply prefix length limit if specified
                    let effective_seq = if prefix_length > 0 && seq1_len > prefix_length {
                        &record_data1.seq[..prefix_length]
                    } else {
                        &record_data1.seq
                    };

                    // Get minimizer hash values
                    fill_minimizer_hashes(
                        effective_seq,
                        kmer_length,
                        window_size,
                        &mut minimizer_buffer1,
                    );
                }

                // Check for minimizer hits in read 2
                if seq2_len >= kmer_length {
                    minimizer_buffer2.clear();

                    // Apply prefix length limit if specified
                    let effective_seq = if prefix_length > 0 && seq2_len > prefix_length {
                        &record_data2.seq[..prefix_length]
                    } else {
                        &record_data2.seq
                    };

                    // Get minimizer hash values
                    fill_minimizer_hashes(
                        effective_seq,
                        kmer_length,
                        window_size,
                        &mut minimizer_buffer2,
                    );
                }

                // Concat the buffers together
                minimizer_buffer1.append(&mut minimizer_buffer2);

                // check_input_matches_index returns true if the matches >= match_threshold
                // If deplete is true, we want to suppress output of sequences that match
                // The NOT is required due to the logic compared to requirement
                // matches_index && deplete = true --> false
                // !matches_index && deplete = false --> true
                // matches_index && !deplete = false --> true
                // !matches_index && !deplete = true --> false
                let should_output = !(check_input_matches_index(
                    minimizer_hashes,
                    &minimizer_buffer1,
                    match_threshold,
                    server_address,
                )
                .await
                    && deplete);

                (should_output, seq1_len, seq2_len)
            })
            .collect();

        // Process results sequentially to maintain order
        for (i, result) in batch_results.into_iter().enumerate() {
            let (should_output, seq1_len, seq2_len) = result.await;
            let record_data1 = &batch1[i];
            let record_data2 = &batch2[i];

            *total_seqs += 2;
            *total_bp += (seq1_len + seq2_len) as u64;

            if should_output {
                // Track output base pairs
                *output_bp += (seq1_len + seq2_len) as u64;

                // Increment output sequence counter (twice, once for each read)
                *output_seq_counter += 2;

                // Format s1 as FASTX to byte buffer and write to appropriate writer
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    &record_data1.id,
                    &record_data1.seq,
                    record_data1.qual.as_deref(),
                    record_data1.format,
                    &mut output_record_buffer,
                    rename,
                    *output_seq_counter - 1,
                );

                if let Some(ref mut w2) = writer2 {
                    // Write read 1 to primary writer
                    writer.write_all(&output_record_buffer)?;

                    // Format s2 as FASTX to byte buffer and write to second writer
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        &record_data2.id,
                        &record_data2.seq,
                        record_data2.qual.as_deref(),
                        record_data2.format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    w2.write_all(&output_record_buffer)?;
                } else {
                    // Interleaved output
                    writer.write_all(&output_record_buffer)?;

                    // Format s2 as FASTX to byte buffer
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        &record_data2.id,
                        &record_data2.seq,
                        record_data2.qual.as_deref(),
                        record_data2.format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    writer.write_all(&output_record_buffer)?;
                }
            } else {
                *filtered_seqs += 2; // Both seqs filtered out
                *filtered_bp += (seq1_len + seq2_len) as u64; // Track filtered base pairs
            }
        }

        // Update spinner and flush periodically
        let elapsed = filtering_start_time.elapsed();
        let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
        let mbp_per_sec = bp_per_sec / 1_000_000.0;

        // Calculate output proportion directly
        let output_seqs = *total_seqs - *filtered_seqs;
        let output_proportion = if *total_seqs > 0 {
            output_seqs as f64 / *total_seqs as f64
        } else {
            0.0
        };

        // Calculate output base pair proportion
        let output_bp_proportion = if *total_bp > 0 {
            *output_bp as f64 / *total_bp as f64
        } else {
            0.0
        };

        // Update spinner message with detailed stats
        spinner.set_message(format!(
            "Retained {}/{} sequences ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            seqs_per_sec,
            mbp_per_sec
        ));

        // Flush writer periodically
        writer.flush()?;
        if let Some(ref mut w2) = writer2 {
            w2.flush()?;
        }

        // Check if we've reached the end of the files
        if reached_end {
            break;
        }
    }

    Ok(())
}

async fn process_interleaved_paired_seqs(
    minimizer_hashes: &Option<FxHashSet<u64>>,
    writer: &mut Box<dyn FastxWriter>,
    mut writer2: Option<&mut Box<dyn FastxWriter>>,
    match_threshold: &MatchThreshold,
    prefix_length: usize,
    kmer_length: usize,
    window_size: usize,
    deplete: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    filtering_start_time: Instant,
    server_address: &Option<String>,
) -> Result<()> {
    // Parse FASTX from stdin
    let mut reader = parse_fastx_stdin()?;
    let mut output_record_buffer = Vec::with_capacity(1024);
    let mut record_counter = 0;

    // Process in batches
    let batch_size = 10000;

    loop {
        // Collect a batch of read pairs with owned data
        let mut batch_pairs = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        // Fill the batch with interleaved pairs
        for _ in 0..batch_size {
            // Read the first record of the pair
            let (record1_id, record1_seq, record1_qual, record1_format) = match reader.next() {
                Some(result) => {
                    record_counter += 1;
                    let record = result?;
                    // Extract all data we need from the record
                    let id = record.id().to_vec();
                    let seq = record.seq().to_vec();
                    let qual = record.qual().map(|q| q.to_vec());
                    let format = record.format();
                    (id, seq, qual, format)
                }
                None => {
                    reached_end = true;
                    break; // End of input
                }
            };

            // Read the second record of the pair
            let (record2_id, record2_seq, record2_qual, record2_format) = match reader.next() {
                Some(result) => {
                    record_counter += 1;
                    let record = result?;
                    let id = record.id().to_vec();
                    let seq = record.seq().to_vec();
                    let qual = record.qual().map(|q| q.to_vec());
                    let format = record.format();
                    (id, seq, qual, format)
                }
                None => {
                    // Check if we have record1 but no record2 (mispaired)
                    return Err(anyhow::anyhow!(
                        "Uneven number of interleaved sequence pairs. Found {} records.",
                        record_counter
                    ));
                }
            };

            // Store the pair in the batch
            batch_pairs.push((
                (record1_id, record1_seq, record1_qual, record1_format),
                (record2_id, record2_seq, record2_qual, record2_format),
            ));
        }

        if batch_pairs.is_empty() {
            break;
        }

        // Process batch in parallel
        let batch_results: Vec<_> = batch_pairs
            .par_iter()
            .map(
                async |(
                    (_record1_id, record1_seq, _record1_qual, _record1_format),
                    (_record2_id, record2_seq, _record2_qual, _record2_format),
                )| {
                    // Pre-allocate buffers for reuse
                    let mut minimizer_buffer1 = Vec::with_capacity(64);
                    let mut minimizer_buffer2 = Vec::with_capacity(64);

                    // Check for minimizer hits in read 1
                    if record1_seq.len() >= kmer_length {
                        // Apply prefix length limit if specified
                        let effective_seq =
                            if prefix_length > 0 && record1_seq.len() > prefix_length {
                                &record1_seq[..prefix_length]
                            } else {
                                &record1_seq
                            };

                        // Get minimizer hash values
                        fill_minimizer_hashes(
                            effective_seq,
                            kmer_length,
                            window_size,
                            &mut minimizer_buffer1,
                        );
                    }

                    // Check for minimizer hits in read 2
                    if record2_seq.len() >= kmer_length {
                        // Apply prefix length limit if specified
                        let effective_seq =
                            if prefix_length > 0 && record2_seq.len() > prefix_length {
                                &record2_seq[..prefix_length]
                            } else {
                                &record2_seq
                            };

                        // Get minimizer hash values
                        fill_minimizer_hashes(
                            effective_seq,
                            kmer_length,
                            window_size,
                            &mut minimizer_buffer2,
                        );
                    }

                    // Concat the buffers together
                    minimizer_buffer1.append(&mut minimizer_buffer2);

                    // check_input_matches_index returns true if the matches >= match_threshold
                    // If deplete is true, we want to suppress output of sequences that match
                    // The NOT is required due to the logic compared to requirement
                    // matches_index && deplete = true --> false
                    // !matches_index && deplete = false --> true
                    // matches_index && !deplete = false --> true
                    // !matches_index && !deplete = true --> false
                    let should_output = !(check_input_matches_index(
                        minimizer_hashes,
                        &minimizer_buffer1,
                        match_threshold,
                        server_address,
                    )
                    .await
                        && deplete);

                    (should_output, record1_seq.len(), record2_seq.len())
                },
            )
            .collect();

        // Process results sequentially to maintain order
        for (i, result) in batch_results.into_iter().enumerate() {
            let (should_output, seq1_len, seq2_len) = result.await;
            let (
                (record1_id, record1_seq, record1_qual, record1_format),
                (record2_id, record2_seq, record2_qual, record2_format),
            ) = &batch_pairs[i];

            *total_seqs += 2;
            *total_bp += (seq1_len + seq2_len) as u64;

            if should_output {
                // Track output base pairs
                *output_bp += (seq1_len + seq2_len) as u64;

                // Increment output sequence counter (twice, once for each seq)
                *output_seq_counter += 2;

                // Format and write record 1
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    record1_id,
                    record1_seq,
                    record1_qual.as_deref(),
                    *record1_format,
                    &mut output_record_buffer,
                    rename,
                    *output_seq_counter - 1,
                );

                if let Some(ref mut w2) = writer2 {
                    // Write read 1 to primary writer
                    writer.write_all(&output_record_buffer)?;

                    // Format and write record 2 to second writer
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        record2_id,
                        record2_seq,
                        record2_qual.as_deref(),
                        *record2_format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    w2.write_all(&output_record_buffer)?;
                } else {
                    // Interleaved output (existing behavior)
                    writer.write_all(&output_record_buffer)?;

                    // Format and write record 2
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        record2_id,
                        record2_seq,
                        record2_qual.as_deref(),
                        *record2_format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    writer.write_all(&output_record_buffer)?;
                }
            } else {
                *filtered_seqs += 2; // Both seqs filtered out
                *filtered_bp += (seq1_len + seq2_len) as u64; // Track filtered base pairs
            }
        }

        // Update spinner and flush periodically
        let elapsed = filtering_start_time.elapsed();
        let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
        let mbp_per_sec = bp_per_sec / 1_000_000.0;

        // Calculate output proportion directly
        let output_seqs = *total_seqs - *filtered_seqs;
        let output_proportion = if *total_seqs > 0 {
            output_seqs as f64 / *total_seqs as f64
        } else {
            0.0
        };

        // Calculate output base pair proportion
        let output_bp_proportion = if *total_bp > 0 {
            *output_bp as f64 / *total_bp as f64
        } else {
            0.0
        };

        // Update spinner message with detailed stats
        spinner.set_message(format!(
            "Retained {}/{} seqs ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            seqs_per_sec,
            mbp_per_sec
        ));

        // Flush writer periodically
        writer.flush()?;
        if let Some(ref mut w2) = writer2 {
            w2.flush()?;
        }

        // Check if we've reached the end of input
        if reached_end {
            break;
        }
    }

    Ok(())
}

/// Push FASTA or FASTQ record to output buffer from component parts
/// Workaround for borrowing misery with interleaved pairs from stdin
fn output_fastx_record_from_parts(
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
    format: Format,
    buffer: &mut Vec<u8>,
    rename: bool,
    seq_number: u64,
) {
    match format {
        Format::Fasta => {
            buffer.push(b'>');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(id);
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(seq);
            buffer.push(b'\n');
        }
        Format::Fastq => {
            buffer.push(b'@');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(id);
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(seq);
            buffer.extend_from_slice(b"\n+\n");
            if let Some(qual_data) = qual {
                buffer.extend_from_slice(qual_data);
            }
            buffer.push(b'\n');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::IndexHeader;
    use crate::index::write_minimizers;
    use tempfile::TempDir;

    #[allow(dead_code)] // Suppress unused warnings
    fn create_test_index() -> (PathBuf, IndexHeader, TempDir) {
        // Create a temporary directory
        let temp_dir = TempDir::new().unwrap();
        let index_path = temp_dir.path().join("test.idx");

        // Create dummy minimizers
        let minimizers: FxHashSet<u64> = [1, 2, 3, 4, 5].iter().cloned().collect();
        let header = IndexHeader::new(5, 3);

        write_minimizers(&minimizers, &header, Some(&index_path)).unwrap();

        // Return the TempDir along with the other values to keep it in scope
        (index_path, header, temp_dir)
    }

    #[test]
    fn test_filter_summary() {
        // Create a sample summary
        let summary = FilterSummary {
            version: "deacon 0.1.0".to_string(),
            index: "test.idx".to_string(),
            input: "test.fastq".to_string(),
            input2: Some("test2.fastq".to_string()),
            output: "output.fastq".to_string(),
            output2: Some("output2.fastq".to_string()),
            k: 31,
            w: 21,
            match_threshold: "1".to_string(),
            prefix_length: 0,
            deplete: false,
            rename: false,
            seqs_in: 100,
            seqs_out: 90,
            seqs_out_proportion: 0.9,
            seqs_removed: 10,
            seqs_removed_proportion: 0.1,
            bp_in: 10000,
            bp_out: 9000,
            bp_out_proportion: 0.9,
            bp_removed: 1000,
            bp_removed_proportion: 0.1,
            time: 1.5,
            seqs_per_second: 66,
            bp_per_second: 6666,
        };

        // Test JSON ser+de
        let json = serde_json::to_string(&summary).unwrap();
        let parsed: FilterSummary = serde_json::from_str(&json).unwrap();

        // Check values
        assert_eq!(parsed.version, "deacon 0.1.0");
        assert_eq!(parsed.seqs_in, 100);
        assert_eq!(parsed.seqs_removed_proportion, 0.1);
        assert_eq!(parsed.seqs_out_proportion, 0.9);
        assert_eq!(parsed.bp_out_proportion, 0.9);
        assert_eq!(parsed.input, "test.fastq");
        assert_eq!(parsed.input2, Some("test2.fastq".to_string()));
        assert_eq!(parsed.output, "output.fastq");
        assert_eq!(parsed.output2, Some("output2.fastq".to_string()));
    }
}
