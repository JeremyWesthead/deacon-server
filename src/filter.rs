use crate::MinimizerSet;
use crate::index::load_minimizers_cached;
use crate::minimizers::KmerHasher;
use crate::minimizers::decode_u64;
use crate::minimizers::decode_u128;
// use crate::index::load_minimizer_hashes;
// use crate::minimizers::fill_minimizer_hashes;
#[cfg(feature = "server")]
use crate::server_common::{FilterRequest, FilterResponse};
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use liblzma::write::XzEncoder;
use needletail::parse_fastx_file;
use needletail::parse_fastx_stdin;
use needletail::parser::Format;
use packed_seq::SeqVec;
use packed_seq::{PackedNSeqVec, u32x8};
use rayon::prelude::*;
#[cfg(feature = "server")]
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use zstd::stream::write::Encoder as ZstdEncoder;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer

#[derive(Clone, Default, Debug)]
pub(crate) struct ProcessingStats {
    pub total_seqs: u64,
    filtered_seqs: u64,
    pub total_bp: u64,
    output_bp: u64,
    filtered_bp: u64,
    output_seq_counter: u64,
    pub last_reported: u64,
}

#[derive(Clone)]
pub(crate) struct Buffers {
    pub packed_nseq: PackedNSeqVec,
    pub positions: Vec<u32>,
    pub minimizers: crate::MinimizerVec,
    pub cache: (simd_minimizers::Cache, Vec<u32x8>, Vec<u32x8>),
}

impl Buffers {
    pub fn new_u64() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            minimizers: crate::MinimizerVec::U64(Vec::new()),
            cache: Default::default(),
        }
    }

    pub fn new_u128() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            minimizers: crate::MinimizerVec::U128(Vec::new()),
            cache: Default::default(),
        }
    }
}

/// Data structure to hold a fastq record
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
                if !(1..=9).contains(&level) {
                    Err(anyhow::anyhow!(
                        "Invalid gzip compression level {}. Must be between 1 and 9.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
            Self::Zstd => {
                if !(1..=22).contains(&level) {
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
            .context(format!("Failed to create output file: {output_path}"))?;

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
    k: u8,
    w: u8,
    abs_threshold: usize,
    rel_threshold: f64,
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

/// Calculate required hits based on absolute and relative thresholds
fn calculate_required_hits(
    total_minimizers: usize,
    abs_threshold: usize,
    rel_threshold: f64,
) -> usize {
    let abs_required = abs_threshold;
    let rel_required = if total_minimizers == 0 {
        0
    } else {
        ((rel_threshold * total_minimizers as f64).round() as usize).max(1)
    };
    abs_required.max(rel_required)
}

/// Check if sequence meets filtering criteria
fn meets_filtering_criteria(
    hit_count: usize,
    total_minimizers: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
) -> bool {
    let required = calculate_required_hits(total_minimizers, abs_threshold, rel_threshold);
    if deplete {
        hit_count < required
    } else {
        hit_count >= required
    }
}

fn get_minimizer_positions_and_values<'s>(
    seq: &'s [u8],
    mut buffers: Buffers,
    prefix_length: usize,
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
) -> Buffers {
    // Apply prefix length limit if specified.
    let seq = if prefix_length > 0 && seq.len() > prefix_length {
        &seq[..prefix_length]
    } else {
        seq
    };
    // Drop trailing newline if present.
    let seq = seq.strip_suffix(b"\n").unwrap_or(seq);

    let Buffers {
        packed_nseq,
        positions,
        minimizers,
        cache,
    } = &mut buffers;

    packed_nseq.seq.clear();
    packed_nseq.ambiguous.clear();
    minimizers.clear();
    positions.clear();

    // Pack the sequence into 2-bit representation.
    packed_nseq.seq.push_ascii(seq);
    packed_nseq.ambiguous.push_ascii(seq);

    // let mut positions = Vec::new();
    let k = kmer_length as usize;
    let w = window_size as usize;
    let m = simd_minimizers::canonical_minimizers(k, w)
        .hasher(&hasher)
        .run_skip_ambiguous_windows_with_buf(packed_nseq.as_slice(), positions, cache);

    // Store k-mer values directly based on variant
    match minimizers {
        crate::MinimizerVec::U64(vec) => {
            vec.extend(m.pos_and_values_u64().map(|(_pos, val)| val));
        }
        crate::MinimizerVec::U128(vec) => {
            vec.extend(m.pos_and_values_u128().map(|(_pos, val)| val));
        }
    }

    buffers
}

fn should_keep_sequence(
    ref_minimizers: &MinimizerSet,
    seq: &[u8],
    kmer_length: u8,
    window_size: u8,
    prefix_length: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
    debug: bool,
) -> (bool, usize, usize, Vec<String>) {
    if seq.len() < kmer_length as usize {
        return (deplete, 0, 0, Vec::new()); // If too short, keep if in deplete mode
    }

    let hasher = KmerHasher::new(kmer_length as usize);
    let buffers = if kmer_length <= 32 {
        Buffers::new_u64()
    } else {
        Buffers::new_u128()
    };

    let buffers = get_minimizer_positions_and_values(
        seq,
        buffers,
        prefix_length,
        kmer_length,
        window_size,
        hasher,
    );
    let minimizers = &buffers.minimizers;
    let num_minimizers = minimizers.len();

    // Count distinct minimizer hits based on variant
    let (hit_count, hit_kmers) = match (minimizers, ref_minimizers) {
        (crate::MinimizerVec::U64(vec), MinimizerSet::U64(set)) => {
            let mut seen_hits = crate::RapidHashSet::default();
            let mut hit_kmers = Vec::new();
            for &minimizer in vec {
                if set.contains(&minimizer) && seen_hits.insert(minimizer) {
                    if debug {
                        let kmer = decode_u64(minimizer, kmer_length);
                        hit_kmers.push(String::from_utf8_lossy(&kmer).to_string());
                    }
                }
            }
            (seen_hits.len(), hit_kmers)
        }
        (crate::MinimizerVec::U128(vec), MinimizerSet::U128(set)) => {
            let mut seen_hits = crate::RapidHashSet::default();
            let mut hit_kmers = Vec::new();
            for &minimizer in vec {
                if set.contains(&minimizer) && seen_hits.insert(minimizer) {
                    if debug {
                        let kmer = decode_u128(minimizer, kmer_length);
                        hit_kmers.push(String::from_utf8_lossy(&kmer).to_string());
                    }
                }
            }
            (seen_hits.len(), hit_kmers)
        }
        _ => panic!("Mismatch between MinimizerVec and MinimizerSet types"),
    };

    (
        meets_filtering_criteria(
            hit_count,
            num_minimizers,
            abs_threshold,
            rel_threshold,
            deplete,
        ),
        hit_count,
        num_minimizers,
        hit_kmers,
    )
}

fn should_keep_pair(
    ref_minimizers: &MinimizerSet,
    seq1: &[u8],
    seq2: &[u8],
    kmer_length: u8,
    window_size: u8,
    prefix_length: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
    debug: bool,
) -> (bool, usize, usize, Vec<String>) {
    // Process both sequences and count distinct hits
    let (_, hit_count1, num_minimizers1, hit_kmers1) = should_keep_sequence(
        ref_minimizers,
        seq1,
        kmer_length,
        window_size,
        prefix_length,
        abs_threshold,
        rel_threshold,
        deplete,
        debug,
    );
    let (_, hit_count2, num_minimizers2, hit_kmers2) = should_keep_sequence(
        ref_minimizers,
        seq2,
        kmer_length,
        window_size,
        prefix_length,
        abs_threshold,
        rel_threshold,
        deplete,
        debug,
    );

    let hit_count = hit_count1 + hit_count2;
    let num_minimizers = num_minimizers1 + num_minimizers2;
    let hit_kmers = [hit_kmers1, hit_kmers2].concat();

    (
        meets_filtering_criteria(
            hit_count,
            num_minimizers,
            abs_threshold,
            rel_threshold,
            deplete,
        ),
        hit_count,
        num_minimizers,
        hit_kmers,
    )
}

/// Given a set of index minimizers and a vector of input minimizers,
/// return a vector of booleans indicating whether each input should be output
pub fn inputs_should_be_output(
    index_minimizers: &MinimizerSet,
    seqs: &Vec<&[u8]>,
    kmer_length: u8,
    window_size: u8,
    prefix_length: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
    debug: bool,
) -> Vec<(bool, Vec<String>, usize)> {
    seqs.par_iter()
        .map(|seq| {
            let (keep, _, _, hit_kmers) = should_keep_sequence(
                index_minimizers,
                seq,
                kmer_length,
                window_size,
                prefix_length,
                abs_threshold,
                rel_threshold,
                deplete,
                debug,
            );

            (keep, hit_kmers, seq.len())
        })
        .collect()
}

/// Given a set of index minimizers and a vector of input minimizers,
/// return a vector of booleans indicating whether each input should be output
pub fn paired_inputs_should_be_output(
    index_minimizers: &MinimizerSet,
    seqs: &Vec<(&[u8], &[u8])>,
    kmer_length: u8,
    window_size: u8,
    prefix_length: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
    debug: bool,
) -> Vec<(bool, Vec<String>, usize, usize)> {
    seqs.par_iter()
        .map(|seq| {
            let (keep, _, _, hit_kmers) = should_keep_pair(
                index_minimizers,
                seq.0,
                seq.1,
                kmer_length,
                window_size,
                prefix_length,
                abs_threshold,
                rel_threshold,
                deplete,
                debug,
            );

            (keep, hit_kmers, seq.0.len(), seq.1.len())
        })
        .collect()
}

// /// Send minimizers to server for checking against index.
// /// Equivalent functionality to `inputs_should_be_output`, but remote
// #[cfg(feature = "server")]
// fn send_all_minimizers_to_server(
//     input_minimizers: Vec<Vec<u64>>,
//     server_address: &str,
//     matches_threshold: &MatchThreshold,
//     deplete: bool,
// ) -> Result<Vec<bool>> {
//     // Create a client to send the minimizers to the server
//     let client = Client::new();

//     // Send the minimizers as a POST request
//     let response = client
//         .post(server_address.to_owned() + "/should_output")
//         .json(&FilterRequest {
//             input: input_minimizers,
//             match_threshold: *matches_threshold,
//             deplete,
//         })
//         .send()?;

//     // Check if the response indicates a match
//     if response.status().is_success() {
//         Ok(response.json::<FilterResponse>()?.should_output)
//     } else {
//         Err(anyhow::anyhow!(
//             "Server returned an error: {}",
//             response.status()
//         ))
//     }
// }

/// Given a set of input minimizers, check if they should be output
/// If index minimizers are provided, check locally.
/// If not, send to server for checking. Requires the `server` feature to be enabled.
pub fn check_inputs_should_be_output(
    index_minimizers: &Option<MinimizerSet>,
    seqs: &Vec<&[u8]>,
    kmer_length: u8,
    window_size: u8,
    prefix_length: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    _server_address: &Option<String>,
    deplete: bool,
    debug: bool,
) -> Vec<(bool, Vec<String>, usize)> {
    // If index minimizers are provided, check if input matches locally
    if let Some(index_minimizers) = index_minimizers {
        inputs_should_be_output(
            index_minimizers,
            seqs,
            kmer_length,
            window_size,
            prefix_length,
            abs_threshold,
            rel_threshold,
            deplete,
            debug,
        )
    } else {
        // Else, send the input minimizers to the server for checking
        // #[cfg(feature = "server")]
        // {
        //     if _server_address.is_none() {
        //         panic!("Server address is required when using the server feature.");
        //     }
        //     let server_address = _server_address.as_ref().map(String::as_str).unwrap();
        //     send_all_minimizers_to_server(
        //         input_minimizers.to_vec(),
        //         server_address,
        //         match_threshold,
        //         deplete,
        //     )
        //     .unwrap_or_else(|e| {
        //         panic!("Error checking input against index: {e}");
        //     })
        // }
        // #[cfg(not(feature = "server"))]
        // {
        //     panic!("Server feature is not enabled. Cannot check input against index.");
        // }
        Vec::new()
    }
}

/// Given a set of input minimizers, check if they should be output
/// If index minimizers are provided, check locally.
/// If not, send to server for checking. Requires the `server` feature to be enabled.
pub fn check_paired_inputs_should_be_output(
    index_minimizers: &Option<MinimizerSet>,
    seqs: &Vec<(&[u8], &[u8])>,
    kmer_length: u8,
    window_size: u8,
    prefix_length: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    _server_address: &Option<String>,
    deplete: bool,
    debug: bool,
) -> Vec<(bool, Vec<String>, usize, usize)> {
    // If index minimizers are provided, check if input matches locally
    if let Some(index_minimizers) = index_minimizers {
        paired_inputs_should_be_output(
            index_minimizers,
            seqs,
            kmer_length,
            window_size,
            prefix_length,
            abs_threshold,
            rel_threshold,
            deplete,
            debug,
        )
    } else {
        // Else, send the input minimizers to the server for checking
        // #[cfg(feature = "server")]
        // {
        //     if _server_address.is_none() {
        //         panic!("Server address is required when using the server feature.");
        //     }
        //     let server_address = _server_address.as_ref().map(String::as_str).unwrap();
        //     send_all_minimizers_to_server(
        //         input_minimizers.to_vec(),
        //         server_address,
        //         match_threshold,
        //         deplete,
        //     )
        //     .unwrap_or_else(|e| {
        //         panic!("Error checking input against index: {e}");
        //     })
        // }
        // #[cfg(not(feature = "server"))]
        // {
        //     panic!("Server feature is not enabled. Cannot check input against index.");
        // }
        Vec::new()
    }
}

/// Run deacon filter with the provided parameters.
#[allow(clippy::too_many_arguments)]
pub fn run(
    minimizers_path: Option<&Path>,
    input_path: &str,
    input2_path: Option<&str>,
    output_path: &str,
    output2_path: Option<&str>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    summary_path: Option<&PathBuf>,
    deplete: bool,
    rename: bool,
    threads: usize,
    compression_level: u8,
    server_address: Option<String>,
    debug: bool,
) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    let tool_version = format!("deacon {version}");

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
    } else if input2_path.is_some() {
        input_type.push_str("paired");
    } else {
        input_type.push_str("single");
    }
    options.push(format!("abs_threshold={abs_threshold}"));
    options.push(format!("rel_threshold={rel_threshold:.3}"));
    if prefix_length > 0 {
        options.push(format!("prefix_length={prefix_length}"));
    }
    if rename {
        options.push("rename".to_string());
    }
    if threads > 0 {
        options.push(format!("threads={threads}"));
    }

    eprintln!(
        "Deacon v{}; mode: {}; input: {}; options: {}",
        version,
        mode,
        input_type,
        options.join(", ")
    );

    // Load minimizers hashes and parse header
    let (minimizer_hashes, header) = load_minimizers_cached(minimizers_path, &server_address)?;

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let load_time = start_time.elapsed();
    eprintln!("Loaded index (k={kmer_length}, w={window_size}) in {load_time:.2?}");

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
            abs_threshold,
            rel_threshold,
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
            debug,
        )?;
    } else if let Some(input2_path) = input2_path {
        process_paired_seqs(
            &minimizer_hashes,
            input_path,
            input2_path,
            &mut writer,
            writer2.as_mut(),
            abs_threshold,
            rel_threshold,
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
            debug,
        )?;
    } else {
        process_single_seqs(
            &minimizer_hashes,
            input_path,
            &mut writer,
            abs_threshold,
            rel_threshold,
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
            debug,
        )?;
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
        "Completed in {total_time:.2?}. Speed: {seqs_per_sec:.0} seqs/s ({mbp_per_sec:.1} Mbp/s)"
    );

    // Build and write a JSON summary if path provided
    if let Some(summary_file) = summary_path {
        // Get number of sequences passing filter
        let seqs_out = total_seqs - filtered_seqs;

        let summary = FilterSummary {
            version: tool_version,
            index: get_summary_index(&minimizers_path, &server_address),
            input: input_path.to_string(),
            input2: input2_path.map(|s| s.to_string()),
            output: output_path.to_string(),
            output2: output2_path.map(|s| s.to_string()),
            k: kmer_length,
            w: window_size,
            abs_threshold,
            rel_threshold,
            prefix_length,
            deplete,
            rename,
            seqs_in: total_seqs,
            seqs_out,
            seqs_out_proportion: output_seq_proportion,
            seqs_removed: filtered_seqs,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp,
            bp_out: output_bp,
            bp_out_proportion: output_bp_proportion,
            bp_removed: filtered_bp,
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
        };

        // Write summary file
        let file = File::create(summary_file)
            .context(format!("Failed to create summary: {summary_file:?}"))?;
        let writer = BufWriter::new(file);

        // Serialise and write the summary JSON
        serde_json::to_writer_pretty(writer, &summary).context("Failed to write summary")?;

        eprintln!("Summary saved to {summary_file:?}");
    }

    Ok(())
}

fn get_summary_index(minimizers_path: &Option<&Path>, server_address: &Option<String>) -> String {
    let index = match minimizers_path {
        Some(path) => path.to_string_lossy().to_string(),
        None => match &server_address {
            None => "No index or server specified".to_string(),
            Some(_addr) => {
                #[cfg(feature = "server")]
                {
                    let client = Client::new();
                    let response = client
                        .get(_addr.to_owned() + "/index_version")
                        .send()
                        .unwrap_or_else(|e| {
                            panic!("Failed to contact server at {}: {e}", _addr);
                        });
                    if response.status().is_success() {
                        _addr.to_owned()
                            + ":"
                            + &response.text().unwrap_or_else(|e| {
                                panic!("Failed to parse server response: {e}");
                            })
                    } else {
                        panic!("Server returned error: {}", response.status())
                    }
                }
                #[cfg(not(feature = "server"))]
                {
                    panic!("Server feature not enabled, cannot use server address");
                }
            }
        },
    };
    index
}

/// Filter a single (unpaired) sequence.
#[allow(clippy::too_many_arguments)]
fn process_single_seqs(
    minimizer_hashes: &Option<MinimizerSet>,
    input_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    kmer_length: u8,
    window_size: u8,
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
    debug: bool,
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

        // Check if minimizers match the index
        // Separated from initial par_iter to allow flexibility with local/server processing
        let batch_should_outputs = check_inputs_should_be_output(
            minimizer_hashes,
            &batch
                .iter()
                .map(|record_data| record_data.seq.as_slice())
                .collect(),
            kmer_length,
            window_size,
            prefix_length,
            abs_threshold,
            rel_threshold,
            server_address,
            deplete,
            debug,
        );

        // Process results sequentially to maintain order
        for (i, (should_output, kmer_hits, seq_len)) in batch_should_outputs.iter().enumerate() {
            let record_data = &batch[i];
            *total_seqs += 1;
            *total_bp += *seq_len as u64;

            if debug && kmer_hits.len() > 0 {
                eprintln!(
                    "DEBUG: {} keep={} kmers=[{}]",
                    String::from_utf8_lossy(&record_data.id),
                    should_output,
                    kmer_hits.join(",")
                );
            }

            if *should_output {
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

/// Filter a pair of sequences
#[allow(clippy::too_many_arguments)]
fn process_paired_seqs(
    minimizer_hashes: &Option<MinimizerSet>,
    input1_path: &str,
    input2_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    mut writer2: Option<&mut Box<dyn FastxWriter>>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    kmer_length: u8,
    window_size: u8,
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
    debug: bool,
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

        // Check if minimizers match the index
        let batch_should_outputs = check_paired_inputs_should_be_output(
            minimizer_hashes,
            &batch1
                .iter()
                .zip(batch2.iter())
                .map(|(record_data1, record_data2)| {
                    (record_data1.seq.as_slice(), record_data2.seq.as_slice())
                })
                .collect(),
            kmer_length,
            window_size,
            prefix_length,
            abs_threshold,
            rel_threshold,
            server_address,
            deplete,
            debug,
        );

        // Process results sequentially to maintain order
        for (i, (should_output, kmer_hits, seq1_len, seq2_len)) in
            batch_should_outputs.iter().enumerate()
        {
            let record_data1 = &batch1[i];
            let record_data2 = &batch2[i];

            *total_seqs += 2;
            *total_bp += (seq1_len + seq2_len) as u64;

            if debug && kmer_hits.len() > 0 {
                eprintln!(
                    "DEBUG: {} {} keep={} kmers=[{}]",
                    String::from_utf8_lossy(&record_data1.id),
                    String::from_utf8_lossy(&record_data2.id),
                    should_output,
                    kmer_hits.join(",")
                );
            }

            if *should_output {
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

/// Filter a pair of interleaved sequences
/// Functionally very similar to `process_paired_seqs`, but handles interleaved input
#[allow(clippy::too_many_arguments)]
fn process_interleaved_paired_seqs(
    minimizer_hashes: &Option<MinimizerSet>,
    writer: &mut Box<dyn FastxWriter>,
    mut writer2: Option<&mut Box<dyn FastxWriter>>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    kmer_length: u8,
    window_size: u8,
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
    debug: bool,
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
                RecordData {
                    id: record1_id,
                    seq: record1_seq,
                    qual: record1_qual,
                    format: record1_format,
                },
                RecordData {
                    id: record2_id,
                    seq: record2_seq,
                    qual: record2_qual,
                    format: record2_format,
                },
            ));
        }

        if batch_pairs.is_empty() {
            break;
        }

        // Check if minimizers match the index
        let batch_should_outputs = check_paired_inputs_should_be_output(
            minimizer_hashes,
            &batch_pairs
                .iter()
                .map(|(record_data1, record_data2)| {
                    (record_data1.seq.as_slice(), record_data2.seq.as_slice())
                })
                .collect(),
            kmer_length,
            window_size,
            prefix_length,
            abs_threshold,
            rel_threshold,
            server_address,
            deplete,
            debug,
        );

        // Process results sequentially to maintain order
        for (i, (should_output, kmer_hits, seq1_len, seq2_len)) in
            batch_should_outputs.iter().enumerate()
        {
            // for (i, result) in batch_results.into_iter().enumerate() {
            let (record1, record2) = &batch_pairs[i];

            *total_seqs += 2;
            *total_bp += (seq1_len + seq2_len) as u64;

            if debug && kmer_hits.len() > 0 {
                eprintln!(
                    "DEBUG: {} {} keep={} kmers=[{}]",
                    String::from_utf8_lossy(&record1.id),
                    String::from_utf8_lossy(&record2.id),
                    should_output,
                    kmer_hits.join(",")
                );
            }

            if *should_output {
                // Track output base pairs
                *output_bp += (seq1_len + seq2_len) as u64;

                // Increment output sequence counter (twice, once for each seq)
                *output_seq_counter += 2;

                // Format and write record 1
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    &record1.id,
                    &record1.seq,
                    record1.qual.as_deref(),
                    record1.format,
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
                        &record2.id,
                        &record2.seq,
                        record2.qual.as_deref(),
                        record2.format,
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
                        &record2.id,
                        &record2.seq,
                        record2.qual.as_deref(),
                        record2.format,
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

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::index::IndexHeader;
//     use crate::index::write_minimizers;
//     use tempfile::TempDir;

//     #[allow(dead_code)] // Suppress unused warnings
//     fn create_test_index() -> (PathBuf, IndexHeader, TempDir) {
//         // Create a temporary directory
//         let temp_dir = TempDir::new().unwrap();
//         let index_path = temp_dir.path().join("test.idx");

//         // Create dummy minimizers
//         let minimizers: FxHashSet<u64> = [1, 2, 3, 4, 5].iter().cloned().collect();
//         let header = IndexHeader::new(5, 3);

//         write_minimizers(&minimizers, &header, Some(&index_path)).unwrap();

//         // Return the TempDir along with the other values to keep it in scope
//         (index_path, header, temp_dir)
//     }

//     #[test]
//     fn test_filter_summary() {
//         // Create a sample summary
//         let summary = FilterSummary {
//             version: "deacon 0.1.0".to_string(),
//             index: "test.idx".to_string(),
//             input: "test.fastq".to_string(),
//             input2: Some("test2.fastq".to_string()),
//             output: "output.fastq".to_string(),
//             output2: Some("output2.fastq".to_string()),
//             k: 31,
//             w: 21,
//             match_threshold: "1".to_string(),
//             prefix_length: 0,
//             deplete: false,
//             rename: false,
//             seqs_in: 100,
//             seqs_out: 90,
//             seqs_out_proportion: 0.9,
//             seqs_removed: 10,
//             seqs_removed_proportion: 0.1,
//             bp_in: 10000,
//             bp_out: 9000,
//             bp_out_proportion: 0.9,
//             bp_removed: 1000,
//             bp_removed_proportion: 0.1,
//             time: 1.5,
//             seqs_per_second: 66,
//             bp_per_second: 6666,
//         };

//         // Test JSON ser+de
//         let json = serde_json::to_string(&summary).unwrap();
//         let parsed: FilterSummary = serde_json::from_str(&json).unwrap();

//         // Check values
//         assert_eq!(parsed.version, "deacon 0.1.0");
//         assert_eq!(parsed.seqs_in, 100);
//         assert_eq!(parsed.seqs_removed_proportion, 0.1);
//         assert_eq!(parsed.seqs_out_proportion, 0.9);
//         assert_eq!(parsed.bp_out_proportion, 0.9);
//         assert_eq!(parsed.input, "test.fastq");
//         assert_eq!(parsed.input2, Some("test2.fastq".to_string()));
//         assert_eq!(parsed.output, "output.fastq");
//         assert_eq!(parsed.output2, Some("output2.fastq".to_string()));
//     }
// }
