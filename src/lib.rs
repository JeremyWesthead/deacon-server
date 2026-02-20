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
pub use index::fetch as index_fetch;
pub use index::{
    INDEX_FORMAT_VERSION, IndexHeader, build as index_build, diff as index_diff,
    dump as index_dump, dump_minimizers, info as index_info, intersect as index_intersect,
    load_minimizers, union as index_union,
};
pub use minimizers::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE};//, compute_minimizer_hashes, fill_minimizer_hashes,
// };

use anyhow::Result;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::hash::BuildHasher;
use std::collections::HashSet;




/// BuildHasher using rapidhash with fixed seed for fast init
#[derive(Clone, Default)]
pub struct FixedRapidHasher;

impl BuildHasher for FixedRapidHasher {
    type Hasher = rapidhash::fast::RapidHasher<'static>;

    fn build_hasher(&self) -> Self::Hasher {
        rapidhash::fast::SeedableState::fixed().build_hasher()
    }
}

/// RapidHashSet using rapidhash with fixed seed for fast init
pub type RapidHashSet<T> = HashSet<T, FixedRapidHasher>;

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
#[derive(Clone)]
pub enum MinimizerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl MinimizerSet {
    pub fn len(&self) -> usize {
        match self {
            MinimizerSet::U64(set) => set.len(),
            MinimizerSet::U128(set) => set.len(),
        }
    }

    pub fn is_u64(&self) -> bool {
        matches!(self, MinimizerSet::U64(_))
    }

    /// Extend with another MinimizerSet (union operation)
    pub fn extend(&mut self, other: Self) {
        match (self, other) {
            (MinimizerSet::U64(self_set), MinimizerSet::U64(other_set)) => {
                self_set.extend(other_set);
            }
            (MinimizerSet::U128(self_set), MinimizerSet::U128(other_set)) => {
                self_set.extend(other_set);
            }
            _ => panic!("Cannot extend U64 set with U128 set or vice versa"),
        }
    }

    /// Remove minimizers from another set (diff operation)
    pub fn remove_all(&mut self, other: &Self) {
        match (self, other) {
            (MinimizerSet::U64(self_set), MinimizerSet::U64(other_set)) => {
                for val in other_set {
                    self_set.remove(val);
                }
            }
            (MinimizerSet::U128(self_set), MinimizerSet::U128(other_set)) => {
                for val in other_set {
                    self_set.remove(val);
                }
            }
            _ => panic!("Cannot remove U128 minimizers from U64 set or vice versa"),
        }
    }

    /// Keep only minimizers present in another set (intersection operation)
    pub fn intersect(&mut self, other: &Self) {
        match (self, other) {
            (MinimizerSet::U64(self_set), MinimizerSet::U64(other_set)) => {
                self_set.retain(|val| other_set.contains(val));
            }
            (MinimizerSet::U128(self_set), MinimizerSet::U128(other_set)) => {
                self_set.retain(|val| other_set.contains(val));
            }
            _ => panic!("Cannot intersect U64 set with U128 set or vice versa"),
        }
    }
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
#[derive(Clone)]
pub enum MinimizerVec {
    U64(Vec<u64>),
    U128(Vec<u128>),
}

impl MinimizerVec {
    pub fn clear(&mut self) {
        match self {
            MinimizerVec::U64(v) => v.clear(),
            MinimizerVec::U128(v) => v.clear(),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            MinimizerVec::U64(v) => v.len(),
            MinimizerVec::U128(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            MinimizerVec::U64(v) => v.is_empty(),
            MinimizerVec::U128(v) => v.is_empty(),
        }
    }
}

/// Match threshold for filtering sequences.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum MatchThreshold {
    Absolute(usize),
    Relative(f64),
}

impl FromStr for MatchThreshold {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(val) = s.parse::<usize>() {
            Ok(MatchThreshold::Absolute(val))
        } else if let Ok(val) = s.parse::<f64>() {
            if val.is_nan() || val.is_sign_negative() || val > 1.0 {
                Err(format!("Relative threshold must be in [0, 1], got: {val}"))
            } else {
                Ok(MatchThreshold::Relative(val))
            }
        } else {
            Err(format!(
                "Invalid threshold format: '{s}'. Expected an integer or a float between [0, 1]"
            ))
        }
    }
}

impl std::fmt::Display for MatchThreshold {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MatchThreshold::Absolute(n) => write!(f, "{n}"),
            MatchThreshold::Relative(p) => write!(f, "{p}"),
        }
    }
}

pub struct FilterConfig {
    /// Minimizer index file path
    pub minimizers_path: PathBuf,

    /// Path to input fastx file (or - for stdin)
    pub input_path: String,

    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<String>,

    /// Path to output fastx file (or - for stdout; detects .gz and .zst)
    pub output_path: String,

    /// Path to optional second output fastx file for paired reads (detects .gz and .zst)
    pub output2_path: Option<String>,

    /// Match threshold for filtering sequences
    pub match_threshold: MatchThreshold,

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
}

impl FilterConfig {
    pub fn new(minimizers_path: PathBuf) -> Self {
        Self {
            minimizers_path,
            input_path: "-".to_string(),
            input2_path: None,
            output_path: "-".to_string(),
            output2_path: None,
            match_threshold: MatchThreshold::Absolute(2),
            prefix_length: 0,
            summary_path: None,
            deplete: false,
            rename: false,
            threads: 0,           // Use all available threads by default
            compression_level: 2, // Default compression level
        }
    }

    pub fn with_input(mut self, input_path: String) -> Self {
        self.input_path = input_path;
        self
    }

    pub fn with_input2(mut self, input2_path: String) -> Self {
        self.input2_path = Some(input2_path);
        self
    }

    pub fn with_output(mut self, output_path: String) -> Self {
        self.output_path = output_path;
        self
    }

    pub fn with_output2(mut self, output2_path: String) -> Self {
        self.output2_path = Some(output2_path);
        self
    }

    pub fn with_match_threshold(mut self, match_threshold: MatchThreshold) -> Self {
        self.match_threshold = match_threshold;
        self
    }

    pub fn with_prefix_length(mut self, prefix_length: usize) -> Self {
        self.prefix_length = prefix_length;
        self
    }

    pub fn with_summary(mut self, summary_path: PathBuf) -> Self {
        self.summary_path = Some(summary_path);
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

    /// Filter with this configuration
    pub fn execute(&self) -> Result<()> {
        filter::run(
            Some(&self.minimizers_path),
            &self.input_path,
            self.input2_path.as_deref(),
            &self.output_path,
            self.output2_path.as_deref(),
            &self.match_threshold,
            self.prefix_length,
            self.summary_path.as_ref(),
            self.deplete,
            self.rename,
            self.threads,
            self.compression_level,
            None,
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

    /// Number of execution threads (0 = auto)
    pub threads: u16,

    /// Suppress per-sequence progress output
    pub quiet: bool,

    /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
    pub entropy_threshold: f32,
}

impl IndexConfig {
    /// Create a new index configuration with the specified input path
    pub fn new(input_path: PathBuf) -> Self {
        Self {
            input_path: input_path,
            kmer_length: DEFAULT_KMER_LENGTH,
            window_size: DEFAULT_WINDOW_SIZE,
            output_path: None,
            threads: 8,
            quiet: false,
            entropy_threshold: 0.0,
        }
    }

    /// Validate k-mer and window size constraints
    pub fn validate(&self) -> Result<()> {
        let k = self.kmer_length as usize;
        let w = self.window_size as usize;

        // Check constraints: k <= 61, k+w <= 96, k+w even (ensures k odd and k+w-1 odd)
        if k > 61 || k + w > 96 || (k + w) % 2 != 0 {
            return Err(anyhow::anyhow!(
                "Invalid k-w combination: k={}, w={}, k+w={} (constraints: k<=61, k+w<=96, k+w even)",
                k,
                w,
                k + w
            ));
        }

        Ok(())
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
    pub fn with_output(mut self, output_path: PathBuf) -> Self {
        self.output_path = Some(output_path);
        self
    }

    /// Set threads
    pub fn with_threads(mut self, threads: u16) -> Self {
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
        index::build(self)
    }
}

// pub fn write_minimizers(
//     minimizers: &FxHashSet<u64>,
//     header: &index::IndexHeader,
//     output_path: Option<&PathBuf>,
// ) -> Result<()> {
//     index::write_minimizers(minimizers, header, output_path)
// }
