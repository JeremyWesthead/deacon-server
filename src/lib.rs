//! # Deacon
//!
//! A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format,
//! built for efficient host depletion (*deacon*-tamination).
//!
//! This crate provides both a library and a binary for filtering nucleotide sequences.
//!
#![doc = include_str!("../README.md")]

// Re-export public functionality
pub mod filter;
pub mod index;
pub mod minimizers;


#[cfg(feature = "server")]
pub mod server;
#[cfg(feature = "server")]
pub mod server_common;

// Re-export the important structures and functions for library users
pub use filter::{FilterSummary, run as run_filter};
pub use index::{
    IndexHeader, build as build_index, diff as diff_index, info as index_info, union as union_index,
};
pub use minimizers::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, compute_minimizer_hashes, fill_minimizer_hashes,
};

use anyhow::Result;
use rustc_hash::FxHashSet;
use std::path::{Path, PathBuf};

pub struct FilterConfig {
    /// Minimizer index file path
    pub minimizers_path: Option<PathBuf>,

    /// Path to input fastx file (or - for stdin)
    pub input_path: String,

    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<String>,

    /// Path to output fastx file (or - for stdout; detects .gz and .zst)
    pub output_path: String,

    /// Path to optional second output fastx file for paired reads (detects .gz and .zst)
    pub output2_path: Option<String>,

    /// Absolute threshold for filtering sequences
    pub abs_threshold: usize,

    /// Relative threshold for filtering sequences (0.0-1.0)
    pub rel_threshold: f64,

    /// Consider only the first N nucleotides per sequence (0 = entire sequence)
    pub prefix_length: usize,

    /// Path to JSON summary file
    pub summary_path: Option<PathBuf>,

    /// Deplete mode (remove sequences WITH matches, original deacon behavior)
    pub deplete: bool,

    /// Replace sequence headers with sequential numbers (1, 2, 3...)
    pub rename: bool,

    /// Number of execution threads (0 = auto)
    pub threads: usize,

    /// Compression level for output files (1-22 for zst, 1-9 for gz)
    pub compression_level: u8,

    /// Debug mode: output sequences with minimizer hits to stderr
    pub debug: bool,

    /// Suppress progress reporting
    pub quiet: bool,

    /// Server address for remote filtering (if using server feature)
    pub server_address: Option<String>,
}

impl FilterConfig {
    pub fn new(minimizers_path: Option<PathBuf>) -> Self {
        Self {
            minimizers_path: minimizers_path,
            input_path: "-".to_string(),
            input2_path: None,
            output_path: "-".to_string(),
            output2_path: None,
            abs_threshold: 2,
            rel_threshold: 0.01,
            prefix_length: 0,
            summary_path: None,
            deplete: false,
            rename: false,
            threads: 0,           // Use all available threads by default
            compression_level: 2, // Default compression level
            debug: false,
            quiet: false,
            server_address: None, // No server address by default
        }
    }

    pub fn with_input<S: Into<String>>(mut self, input_path: S) -> Self {
        self.input_path = input_path.into();
        self
    }

    pub fn with_input2<S: Into<String>>(mut self, input2_path: S) -> Self {
        self.input2_path = Some(input2_path.into());
        self
    }

    pub fn with_output<S: Into<String>>(mut self, output_path: S) -> Self {
        self.output_path = output_path.into();
        self
    }

    pub fn with_output2<S: Into<String>>(mut self, output2_path: S) -> Self {
        self.output2_path = Some(output2_path.into());
        self
    }

    pub fn with_abs_threshold(mut self, abs_threshold: usize) -> Self {
        self.abs_threshold = abs_threshold;
        self
    }

    pub fn with_rel_threshold(mut self, rel_threshold: f64) -> Self {
        self.rel_threshold = rel_threshold;
        self
    }

    pub fn with_prefix_length(mut self, prefix_length: usize) -> Self {
        self.prefix_length = prefix_length;
        self
    }

    pub fn with_summary<P: AsRef<Path>>(mut self, summary_path: P) -> Self {
        self.summary_path = Some(summary_path.as_ref().to_path_buf());
        self
    }

    pub fn with_deplete(mut self, deplete: bool) -> Self {
        self.deplete = deplete;
        self
    }

    pub fn with_rename(mut self, rename: bool) -> Self {
        self.rename = rename;
        self
    }

    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn with_compression_level(mut self, compression_level: u8) -> Self {
        self.compression_level = compression_level;
        self
    }

    pub fn with_debug(mut self, debug: bool) -> Self {
        self.debug = debug;
        self
    }

    pub fn with_quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    pub fn with_server_address(mut self, server_address: String) -> Self {
        self.server_address = Some(server_address);
        self
    }

    /// Filter with this configuration
    pub fn execute(&self) -> Result<()> {
        filter::run(
            self.minimizers_path.as_ref(),
            &self.input_path,
            self.input2_path.as_deref(),
            &self.output_path,
            self.output2_path.as_deref(),
            self.abs_threshold,
            self.rel_threshold,
            self.prefix_length,
            self.summary_path.as_ref(),
            self.deplete,
            self.rename,
            self.threads,
            self.compression_level,
            self.debug,
            self.quiet,
            &self.server_address,
        )
    }
}

pub struct IndexConfig {
    /// Path to input fastx file
    pub input_path: PathBuf,

    /// K-mer length used for indexing
    pub kmer_length: u8,

    /// Minimizer window size used for indexing
    pub window_size: u8,

    /// Path to output file (None for stdout)
    pub output_path: Option<PathBuf>,

    /// Hash table pre-allocation capacity in millions
    pub capacity_millions: usize,

    /// Number of execution threads (0 = auto)
    pub threads: usize,

    /// Suppress per-sequence progress output
    pub quiet: bool,

    /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
    pub entropy_threshold: f32,
}

impl IndexConfig {
    /// Create a new index configuration with the specified input path
    pub fn new<P: AsRef<Path>>(input_path: P) -> Self {
        Self {
            input_path: input_path.as_ref().to_path_buf(),
            kmer_length: DEFAULT_KMER_LENGTH,
            window_size: DEFAULT_WINDOW_SIZE,
            output_path: None,
            capacity_millions: 400,
            threads: 8,
            quiet: false,
            entropy_threshold: 0.0,
        }
    }

    /// Set k-mer length
    pub fn with_kmer_length(mut self, kmer_length: u8) -> Self {
        self.kmer_length = kmer_length;
        self
    }

    /// Set window size
    pub fn with_window_size(mut self, window_size: u8) -> Self {
        self.window_size = window_size;
        self
    }

    /// Set output path
    pub fn with_output<P: AsRef<Path>>(mut self, output_path: P) -> Self {
        self.output_path = Some(output_path.as_ref().to_path_buf());
        self
    }

    /// Set hash table capacity in millions
    pub fn with_capacity_millions(mut self, capacity_millions: usize) -> Self {
        self.capacity_millions = capacity_millions;
        self
    }

    /// Set threads
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set quiet mode
    pub fn with_quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    /// Set threshold for scaled entropy filtering at indexing time
    pub fn with_entropy_threshold(mut self, threshold: f32) -> Self {
        self.entropy_threshold = threshold;
        self
    }

    /// Execute index build with this configuration
    pub fn execute(&self) -> Result<()> {
        build_index(
            &self.input_path,
            self.kmer_length,
            self.window_size,
            self.output_path.clone(),
            self.capacity_millions,
            self.threads,
            self.quiet,
            self.entropy_threshold,
        )
    }
}

pub fn load_minimizers(path: &PathBuf) -> Result<(Option<FxHashSet<u64>>, index::IndexHeader)> {
    index::load_minimizer_hashes(&Some(path), &None)
}

pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &index::IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    index::write_minimizers(minimizers, header, output_path)
}
