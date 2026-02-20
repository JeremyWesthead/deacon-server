use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use deacon::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, IndexConfig, MatchThreshold, index_diff, index_dump,
    index_fetch, index_info, index_intersect, index_union, run_filter,
};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create and compose minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommands,
    },
    /// Keep or discard fastx records with sufficient minimizer hits to the index
    Filter {
        /// Path to minimizer index file
        index: PathBuf,

        /// Optional path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Optional path to second paired fastx file (or - for interleaved stdin)
        input2: Option<String>,

        /// Path to output fastx file (or - for stdout; detects .gz and .zst)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Optional path to second paired output fastx file (detects .gz and .zst)
        #[arg(short = 'O', long = "output2")]
        output2: Option<String>,

        /// Minimum absolute number of minimizer hits for a match
        #[arg(short = 'a', long = "abs-threshold", default_value_t = 2, value_parser = clap::value_parser!(u8).range(1..))]
        abs_threshold: u8,

        /// Minimum relative proportion (0.0-1.0) of minimizer hits for a match
        #[arg(short = 'r', long = "rel-threshold", default_value_t = 0.01)]
        rel_threshold: f64,

        /// Search only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'p', long = "prefix-length", default_value_t = 0)]
        prefix_length: usize,

        /// Discard matching sequences (invert filtering behaviour)
        #[arg(short = 'd', long = "deplete", default_value_t = false)]
        deplete: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'r', long = "rename", default_value_t = false)]
        rename: bool,

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Output compression level (1-9 for gz & xz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 2)]
        compression_level: u8,

        /// Output sequences with minimizer hits to stderr
        #[arg(long = "debug", default_value_t = false)]
        debug: bool,
    },

    Server {
        /// Path to minimizer index file
        index: PathBuf,

        /// Port to run the server on
        #[arg(short = 'p', long = "port", default_value_t = 8888)]
        port: u16,
    },

    /// Alternate version of Filter, swapping local compute for passing to a server
    /// which has the index pre-loaded. Will inevitably be slower than local filtering,
    /// but saves on index loading. Better used for cases of small input + large index
    Client {
        /// Server address to connect to (including port)
        server_address: String,

        /// Optional path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Optional path to second paired fastx file (or - for interleaved stdin)
        input2: Option<String>,

        /// Path to output fastx file (or - for stdout; detects .gz and .zst)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Optional path to second paired output fastx file (detects .gz and .zst)
        #[arg(short = 'O', long = "output2")]
        output2: Option<String>,

        /// Minimum absolute number of minimizer hits for a match
        #[arg(short = 'a', long = "abs-threshold", default_value_t = 2, value_parser = clap::value_parser!(u8).range(1..))]
        abs_threshold: u8,

        /// Minimum relative proportion (0.0-1.0) of minimizer hits for a match
        #[arg(short = 'r', long = "rel-threshold", default_value_t = 0.01)]
        rel_threshold: f64,

        /// Search only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'p', long = "prefix-length", default_value_t = 0)]
        prefix_length: usize,

        /// Discard matching sequences (invert filtering behaviour)
        #[arg(short = 'd', long = "deplete", default_value_t = false)]
        deplete: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'r', long = "rename", default_value_t = false)]
        rename: bool,

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Output compression level (1-9 for gz & xz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 2)]
        compression_level: u8,

        /// Output sequences with minimizer hits to stderr
        #[arg(long = "debug", default_value_t = false)]
        debug: bool,
    },
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Index minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (or - for stdin; supports gz, zst and xz compression)
        input: PathBuf,

        /// K-mer length used for indexing (k+w-1 must be <= 96 and odd)
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: u8,

        /// Minimizer window size used for indexing
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,

        /// Suppress sequence header output
        #[arg(short = 'q', long = "quiet")]
        quiet: bool,

        /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
        #[arg(short = 'e', long = "entropy-threshold", default_value = "0.0")]
        entropy_threshold: f32,
    },
    /// Combine multiple minimizer indexes (A ∪ B…)
    Union {
        /// Path(s) to one or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Intersect multiple minimizer indexes (A ∩ B…)
    Intersect {
        /// Path(s) to two or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Subtract minimizers in one index from another (A - B)
    Diff {
        /// Path to first index file
        #[arg(required = true)]
        first: PathBuf,

        /// Path to second index file or FASTX file (or - for stdin when using FASTX)
        #[arg(required = true)]
        second: PathBuf,

        /// K-mer length (required if second argument is FASTX file, 1-32)
        #[arg(short = 'k', long = "kmer-length", value_parser = clap::value_parser!(u8).range(1..=32))]
        kmer_length: Option<u8>,

        /// Window size (required if second argument is FASTX file)
        #[arg(short = 'w', long = "window-size")]
        window_size: Option<u8>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Dump minimizer index to fasta
    Dump {
        /// Path to index file
        index: PathBuf,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Show index information
    Info {
        /// Path to index file
        index: PathBuf,
    },
    /// Fetch a pre-built index from remote storage
    Fetch {
        /// Index name (e.g., panhuman-1)
        #[arg(default_value = "panhuman-1")]
        index_name: String,

        /// K-mer length
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: u8,

        /// Minimizer window size
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Path to output file (default: ./)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
}
fn main() -> Result<()> {
    // Check we have either AVX2 or NEON
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        eprintln!(
            "Warning: SIMD acceleration is unavailable. For best performance, compile with `cargo build --release -C target-cpu=native`"
        );
    }

    let cli = Cli::parse();

    match &cli.command {
        Commands::Index { command } => match command {
            IndexCommands::Build {
                input,
                kmer_length,
                window_size,
                output,
                threads,
                quiet,
                entropy_threshold,
            } => {
                let config = IndexConfig {
                    input_path: input.clone(),
                    kmer_length: *kmer_length,
                    window_size: *window_size,
                    output_path: output.clone(),
                    threads: *threads,
                    quiet: *quiet,
                    entropy_threshold: *entropy_threshold,
                };
                config
                    .execute()
                    .context("Failed to run index build command")?;
            }
            IndexCommands::Info { index } => {
                index_info(index).context("Failed to run index info command")?;
            }
            IndexCommands::Fetch {
                index_name,
                kmer_length,
                window_size,
                output,
            } => {
                index_fetch(index_name, *kmer_length, *window_size, output.as_deref())
                    .context("Failed to run index fetch command")?;
            }
            IndexCommands::Union { inputs, output } => {
                index_union(inputs, output.as_deref())
                    .context("Failed to run index union command")?;
            }
            IndexCommands::Intersect { inputs, output } => {
                index_intersect(inputs, output.as_deref())
                    .context("Failed to run index intersect command")?;
            }
            IndexCommands::Diff {
                first,
                second,
                kmer_length,
                window_size,
                output,
                threads,
            } => {
                index_diff(
                    first,
                    second,
                    *kmer_length,
                    *window_size,
                    *threads,
                    output.as_deref(),
                )
                .context("Failed to run index diff command")?;
            }
            IndexCommands::Dump { index, output } => {
                index_dump(index, output.as_deref()).context("Failed to run index dump command")?;
            }
        },
        Commands::Filter {
            index: minimizers,
            input,
            input2,
            output,
            output2,
            abs_threshold,
            rel_threshold,
            prefix_length,
            summary,
            deplete,
            rename,
            threads,
            compression_level,
            debug,
        } => {
            // Validate output2 usage
            if output2.is_some() && input2.is_none() {
                eprintln!(
                    "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                );
            }

            run_filter(
                Some(minimizers),
                input,
                input2.as_deref(),
                output,
                output2.as_deref(),
                *abs_threshold as usize,
                *rel_threshold,
                *prefix_length,
                summary.as_ref(),
                *deplete,
                *rename,
                *threads,
                *compression_level,
                None,
                *debug,
            )
            .context("Failed to run filter command")?;
        }
        Commands::Server { index, port } => {
            #[cfg(feature = "server")]
            {
                // Server needs to run async, so spawn an async runtime to run it
                let rt = tokio::runtime::Runtime::new().unwrap();
                rt.block_on(async {
                    eprintln!("Loading server!");
                    deacon::server::run_server(index.clone(), *port).await;
                });
            }
            #[cfg(not(feature = "server"))]
            {
                eprintln!(
                    "Server functionality is not enabled in this build. Please compile with the 'server' feature: `cargo build --features server`"
                );
                // Suppress dead code warning so this compiles without issue when server is not enabled
                let _ = (index, port);
                std::process::exit(1);
            }
        }
        Commands::Client {
            server_address,
            input,
            input2,
            output,
            output2,
            abs_threshold,
            rel_threshold,
            prefix_length,
            summary,
            deplete,
            rename,
            threads,
            compression_level,
            debug,
        } => {
            #[cfg(feature = "server")]
            {
                // Validate output2 usage

                use std::path::Path;
                if output2.is_some() && input2.is_none() {
                    eprintln!(
                        "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                    );
                }
                let dummy_index: Option<&Path> = None;
                run_filter(
                    dummy_index,
                    input,
                    input2.as_deref(),
                    output,
                    output2.as_deref(),
                    *abs_threshold as usize,
                    *rel_threshold,
                    *prefix_length,
                    summary.as_ref(),
                    *deplete,
                    *rename,
                    *threads,
                    *compression_level,
                    Some(server_address.to_string()),
                    *debug,
                )
                .context("Failed to run filter with client functionality")?;
            }
            #[cfg(not(feature = "server"))]
            {
                eprintln!(
                    "Client functionality is not enabled in this build. Please compile with the 'server' feature: `cargo build --features server`"
                );
                // Suppress dead code warning so this compiles without issue when server is not enabled
                let _ = (
                    server_address,
                    input,
                    input2,
                    output,
                    output2,
                    abs_threshold,
                    rel_threshold,
                    prefix_length,
                    summary,
                    deplete,
                    rename,
                    threads,
                    compression_level,
                    debug,
                );
                std::process::exit(1);
            }
        }
    }

    Ok(())
}
