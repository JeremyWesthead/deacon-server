use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use needletail::{parse_fastx_file, parse_fastx_stdin};

use crate::server_common::get_sever_index_header;

/// Serialisable index file
#[derive(Serialize, Deserialize, Debug)]
pub struct Index {
    pub header: IndexHeader,
    pub minimizers: FxHashSet<u64>,
}

/// Serialisable header for the index file
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct IndexHeader {
    pub format_version: u8,
    pub kmer_length: u8,
    pub window_size: u8,
}

impl IndexHeader {
    pub fn new(kmer_length: usize, window_size: usize) -> Self {
        IndexHeader {
            format_version: 2,
            kmer_length: kmer_length as u8,
            window_size: window_size as u8,
        }
    }

    /// Validate header
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.format_version != 2 {
            return Err(anyhow::anyhow!(
                "Unsupported index format version: {}",
                self.format_version
            ));
        }

        Ok(())
    }

    /// Get k
    pub fn kmer_length(&self) -> usize {
        self.kmer_length as usize
    }

    /// Get w
    pub fn window_size(&self) -> usize {
        self.window_size as usize
    }
}

/// Load just the header and count from an index file
pub fn load_header_and_count<P: AsRef<Path>>(path: &P) -> Result<(IndexHeader, usize)> {
    let file =
        File::open(path).context(format!("Failed to open index file {:?}", path.as_ref()))?;
    let mut reader = BufReader::new(file);

    // Deserialise header
    let header: IndexHeader = decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers
    let count: usize = decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialise minimizer count")?;

    Ok((header, count))
}

/// Load the hashes without spiking memory usage with an extra vec
pub async fn load_minimizer_hashes<P: AsRef<Path>>(
    path_option: &Option<P>,
    server_address_option: &Option<String>,
) -> Result<(Option<FxHashSet<u64>>, IndexHeader)> {
    if let Some(path) = path_option {
        let file =
            File::open(path).context(format!("Failed to open index file {:?}", path.as_ref()))?;
        let mut reader = BufReader::new(file);

        // Deserialise header
        let header: IndexHeader = decode_from_std_read(&mut reader, bincode::config::standard())
            .context("Failed to deserialise index header")?;
        header.validate()?;

        // Deserialise the count of minimizers so we can init a FxHashSet with the right capacity
        let count: usize = decode_from_std_read(&mut reader, bincode::config::standard())
            .context("Failed to deserialise minimizer count")?;

        // Pre-allocate FxHashSet with correct capacity
        let mut minimizers = FxHashSet::with_capacity_and_hasher(count, Default::default());

        // Populate FxHashSet
        for _ in 0..count {
            let hash: u64 = decode_from_std_read(&mut reader, bincode::config::standard())
                .context("Failed to deserialise minimizer hash")?;
            minimizers.insert(hash);
        }

        Ok((Some(minimizers), header))
    } else {
        // If no path is provided, check if a server adress was given to populate the header
        #[cfg(feature = "server")]
        {
            if let Some(server) = server_address_option {
                return Ok((None, get_sever_index_header(&server.to_string()).await?));
            } else {
                return Err(anyhow::anyhow!(
                    "No server address provided for running in server mode"
                ));
            }
        }
        #[cfg(not(feature = "server"))]
        {
            return Err(anyhow::anyhow!(
                "Server feature is not enabled. Cannot run without an index."
            ));
        }
    }
}

/// Helper function to write minimizers to output file or stdout
pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    // Create writer based on output path
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        if path.to_string_lossy() == "-" {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            Box::new(BufWriter::new(
                File::create(path).context("Failed to create output file")?,
            ))
        }
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    // Serialise header and minimizers
    let mut writer = BufWriter::new(writer);
    encode_into_std_write(header, &mut writer, bincode::config::standard())
        .context("Failed to serialise index header")?;

    // Serialise the count of minimizers first
    let count = minimizers.len();
    encode_into_std_write(&count, &mut writer, bincode::config::standard())
        .context("Failed to serialise minimizer count")?;

    // Serialise each minimizer directly
    for &hash in minimizers {
        encode_into_std_write(&hash, &mut writer, bincode::config::standard())
            .context("Failed to serialise minimizer hash")?;
    }
    Ok(())
}

/// Build an index of minimizers from a fastx file
pub fn build<P: AsRef<Path>>(
    input: P,
    kmer_length: usize,
    window_size: usize,
    output: Option<PathBuf>,
    capacity_millions: usize,
    threads: usize,
) -> Result<()> {
    let start_time = Instant::now();
    let path = input.as_ref();

    let version: String = env!("CARGO_PKG_VERSION").to_string();

    // Build options string similar to filter
    let mut options = Vec::<String>::new();
    options.push(format!("capacity={}M", capacity_millions));
    if threads > 0 {
        options.push(format!("threads={}", threads));
    }

    eprintln!(
        "Deacon v{}; mode: build; input: single; options: {}",
        version,
        options.join(", ")
    );

    // Ensure l = k + w - 1 is odd so that canonicalisation tie breaks work correctly
    let l = kmer_length + window_size - 1;
    if l % 2 == 0 {
        return Err(anyhow::anyhow!(
            "Constraint violated: k + w - 1 must be odd (k={}, w={})",
            kmer_length,
            window_size
        ));
    }

    // Configure thread pool if specified (non-zero)
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to initialize thread pool")?;
    }

    // Create a fastx reader using needletail (handles gzip automatically)
    let mut reader = parse_fastx_file(path).context("Failed to open input file")?;

    // Init FxHashSet with user-specified capacity
    let capacity = capacity_millions * 1_000_000;
    let mut all_minimizers: FxHashSet<u64> =
        FxHashSet::with_capacity_and_hasher(capacity, Default::default());

    eprintln!("Building index (k={}, w={})", kmer_length, window_size);

    let mut seq_count = 0;
    let mut total_bp = 0;

    // Process sequences in batches for parallelization
    let batch_size = 10000;

    loop {
        // Collect a batch of sequences
        let mut batch = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        for _ in 0..batch_size {
            if let Some(record_result) = reader.next() {
                match record_result {
                    Ok(record) => {
                        let seq_data = record.seq().to_vec();
                        let id = record.id().to_vec();
                        batch.push((seq_data, id));
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
        let batch_results: Vec<_> = batch
            .par_iter()
            .map(|(seq_data, _id)| {
                // Compute minimizer hashes for this sequence
                crate::minimizers::compute_minimizer_hashes(seq_data, kmer_length, window_size)
            })
            .collect();

        // Merge results sequentially (to avoid concurrent modification of hash set)
        for (i, hashes) in batch_results.iter().enumerate() {
            let (seq_data, id) = &batch[i];

            all_minimizers.extend(hashes.iter());

            seq_count += 1;
            total_bp += seq_data.len();

            let id_str = std::str::from_utf8(id).unwrap_or("unknown");
            eprintln!(
                "  {} ({}bp), total minimizers: {}",
                id_str,
                seq_data.len(),
                all_minimizers.len()
            );
        }

        // Check if we've reached the end of the file
        if reached_end {
            break;
        }
    }

    eprintln!(
        "Indexed {} minimizers from {} sequence(s) ({}bp)",
        all_minimizers.len(),
        seq_count,
        total_bp
    );

    let header = IndexHeader::new(kmer_length, window_size);

    // Write to output path or stdout
    write_minimizers(&all_minimizers, &header, output.as_ref())?;

    let total_time = start_time.elapsed();
    eprintln!("Completed in {:.2?}", total_time);

    Ok(())
}

/// Stream minimizers from a FASTX file or stdin and remove those present in first_minimizers
fn stream_diff_fastx<P: AsRef<Path>>(
    fastx_path: P,
    kmer_length: usize,
    window_size: usize,
    first_header: &IndexHeader,
    first_minimizers: &mut FxHashSet<u64>,
) -> Result<(usize, usize)> {
    let path = fastx_path.as_ref();

    // Validate parameters match the first index
    if kmer_length != first_header.kmer_length() || window_size != first_header.window_size() {
        return Err(anyhow::anyhow!(
            "FASTX parameters (k={}, w={}) must match first index (k={}, w={})",
            kmer_length,
            window_size,
            first_header.kmer_length(),
            first_header.window_size()
        ));
    }

    if path.to_string_lossy() == "-" {
        eprintln!(
            "Second index: processing FASTX from stdin (k={}, w={})…",
            kmer_length, window_size
        );
    } else {
        eprintln!(
            "Second index: processing FASTX from file (k={}, w={})…",
            kmer_length, window_size
        );
    }

    // Create a fastx reader
    let mut reader = if path.to_string_lossy() == "-" {
        parse_fastx_stdin().context("Failed to parse FASTX from stdin")?
    } else {
        parse_fastx_file(path).context("Failed to open FASTX file")?
    };

    let mut seq_count = 0;
    let mut total_bp = 0;
    let mut removed_count = 0;
    let mut last_reported_gb = 0;

    // Process sequences in batches for parallelization
    let batch_size = 1000;

    loop {
        // Collect a batch of sequences
        let mut batch = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        for _ in 0..batch_size {
            if let Some(record_result) = reader.next() {
                match record_result {
                    Ok(record) => {
                        let seq_data = record.seq().to_vec();
                        let id = record.id().to_vec();
                        batch.push((seq_data, id));
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
        let batch_results: Vec<_> = batch
            .par_iter()
            .map(|(seq_data, _id)| {
                // Compute minimizer hashes for this sequence
                crate::minimizers::compute_minimizer_hashes(seq_data, kmer_length, window_size)
            })
            .collect();

        // Stream removal: remove minimizers as we find them
        for (i, hashes) in batch_results.iter().enumerate() {
            let (seq_data, _id) = &batch[i];

            // Remove matching minimizers from first_minimizers immediately
            for &hash in hashes {
                if first_minimizers.remove(&hash) {
                    removed_count += 1;
                }
            }

            seq_count += 1;
            total_bp += seq_data.len();

            let current_gb = total_bp / 1_000_000_000;
            if current_gb > last_reported_gb {
                eprintln!(
                    "  Processed {} sequences ({}bp), removed {} minimizers",
                    seq_count, total_bp, removed_count
                );
                last_reported_gb = current_gb;
            }
        }

        // Check if we've reached the end of the file
        if reached_end {
            break;
        }
    }

    eprintln!(
        "Processed {} sequences ({}bp) from FASTX file",
        seq_count, total_bp
    );

    Ok((seq_count, total_bp))
}

/// Compute the set difference between two minimizer indexes (A - B)
pub async fn diff<P: AsRef<Path>>(
    first: P,
    second: P,
    kmer_length: Option<usize>,
    window_size: Option<usize>,
    output: Option<&PathBuf>,
) -> Result<()> {
    let start_time = Instant::now();
    let first = Some(first);
    let second = Some(second);

    // Load first file (always an index)
    let (mut first_minimizers, header) = load_minimizer_hashes(&first, &None).await?;
    if first_minimizers.is_none() {
        return Err(anyhow::anyhow!("Failed to load first index file"));
    }
    let mut first_minimizers = first_minimizers.unwrap();
    eprintln!("First index: loaded {} minimizers", first_minimizers.len());

    // Guess if second file is an index or FASTX file
    let second_minimizers = if let (Some(k), Some(w)) = (kmer_length, window_size) {
        // Second file is a FASTX file - stream diff with provided k, w
        let before_count = first_minimizers.len();
        let (_seq_count, _total_bp) =
            stream_diff_fastx(&second.unwrap(), k, w, &header, &mut first_minimizers)?;

        // Report results
        eprintln!(
            "Removed {} minimizers, {} remaining",
            before_count - first_minimizers.len(),
            first_minimizers.len()
        );

        write_minimizers(&first_minimizers, &header, output)?;

        let total_time = start_time.elapsed();
        eprintln!("Completed difference operation in {:.2?}", total_time);

        return Ok(());
    } else {
        // Try to load as index file first
        if let Ok((second_minimizers, second_header)) = load_minimizer_hashes(&second, &None).await
        {
            if second_minimizers.is_none() {
                return Err(anyhow::anyhow!("Failed to load second index file"));
            }
            let second_minimizers = second_minimizers.unwrap();
            // Second file is an index file
            eprintln!(
                "Second index: loaded {} minimizers",
                second_minimizers.len()
            );

            // Validate compatible headers
            if second_header.kmer_length() != header.kmer_length()
                || second_header.window_size() != header.window_size()
            {
                return Err(anyhow::anyhow!(
                    "Incompatible headers: second index has k={}, w={}, but first index has k={}, w={}",
                    second_header.kmer_length(),
                    second_header.window_size(),
                    header.kmer_length(),
                    header.window_size()
                ));
            }

            second_minimizers
        } else {
            // Second file is not a valid index, treat as FASTX file
            // Use k and w from first index header and do a streaming diff
            let k = header.kmer_length();
            let w = header.window_size();
            // eprintln!(
            //     "Second index is not valid, treating as FASTX and extracting minimizers (k={}, w={})",
            //     k, w
            // );

            // Count minimizers before diff
            let before_count = first_minimizers.len();

            let (_seq_count, _total_bp) =
                stream_diff_fastx(&second.unwrap(), k, w, &header, &mut first_minimizers)?;

            // Report results
            eprintln!(
                "Removed {} minimizers, {} remaining",
                before_count - first_minimizers.len(),
                first_minimizers.len()
            );

            write_minimizers(&first_minimizers, &header, output)?;

            let total_time = start_time.elapsed();
            eprintln!("Completed difference operation in {:.2?}", total_time);

            return Ok(());
        }
    };

    // Handle straightforward index-to-index diffing
    // Count minimizers before diff
    let before_count = first_minimizers.len();

    // Remove all hashes in second_minimizers from first_minimizers
    for hash in &second_minimizers {
        first_minimizers.remove(hash);
    }

    // Report results
    eprintln!(
        "Removed {} minimizers, {} remaining",
        before_count - first_minimizers.len(),
        first_minimizers.len()
    );

    write_minimizers(&first_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    eprintln!("Completed diff operation in {:.2?}", total_time);

    Ok(())
}

/// Show info about an index
pub async fn info<P: AsRef<Path>>(index_path: P) -> Result<()> {
    let start_time = Instant::now();

    // Load index file
    let (minimizers, header) = load_minimizer_hashes(&Some(index_path), &None).await?;
    if minimizers.is_none() {
        return Err(anyhow::anyhow!("Failed to load index file"));
    }
    let minimizers = minimizers.unwrap();
    // Show index info
    eprintln!("Index information:");
    eprintln!("  Format version: {}", header.format_version);
    eprintln!("  K-mer length (k): {}", header.kmer_length());
    eprintln!("  Window size (w): {}", header.window_size());
    eprintln!("  Distinct minimizer count: {}", minimizers.len());

    let total_time = start_time.elapsed();
    eprintln!("Retrieved index info in {:.2?}", total_time);

    Ok(())
}

/// Combine minimizer indexes (set union)
pub async fn union<P: AsRef<Path>>(inputs: &[P], output: Option<&PathBuf>) -> Result<()> {
    let start_time = Instant::now();
    // Check input files
    if inputs.is_empty() {
        return Err(anyhow::anyhow!(
            "No input files provided for union operation"
        ));
    }

    // Read all headers first to determine total capacity needed
    let mut headers_and_counts = Vec::new();
    let mut total_capacity = 0;

    for path in inputs {
        let (header, count) = load_header_and_count(path)?;
        total_capacity += count;
        headers_and_counts.push((header, count));
    }

    // Get header from first file for output
    let header = &headers_and_counts[0].0;

    eprintln!(
        "Performing union of indexes (k={}, w={})",
        header.kmer_length(),
        header.window_size()
    );
    eprintln!(
        "Pre-allocating capacity for {} minimizers across {} indexes",
        total_capacity,
        inputs.len()
    );

    // Verify all headers are compatible
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length() != header.kmer_length()
            || file_header.window_size() != header.window_size()
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={}",
                i,
                file_header.kmer_length(),
                file_header.window_size(),
                header.kmer_length(),
                header.window_size()
            ));
        }
    }

    // Pre-allocate hash set with total capacity to avoid resizing
    let mut all_minimizers: FxHashSet<u64> =
        FxHashSet::with_capacity_and_hasher(total_capacity, Default::default());

    // Now load and merge all indexes
    for (i, path) in inputs.iter().enumerate() {
        let (minimizers, _) = load_minimizer_hashes(&Some(path), &None).await?;
        let before_count = all_minimizers.len();
        if minimizers.is_none() {
            return Err(anyhow::anyhow!(
                "Failed to load index file {:?}",
                path.as_ref()
            ));
        }
        let minimizers = minimizers.unwrap();

        // Merge minimizers (set union)
        all_minimizers.extend(minimizers);

        let expected_count = headers_and_counts[i].1;
        eprintln!(
            "Index {}: expected {} minimizers, added {} new, total: {}",
            i + 1,
            expected_count,
            all_minimizers.len() - before_count,
            all_minimizers.len()
        );
    }

    write_minimizers(&all_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    eprintln!(
        "United {} indexes with {} total minimizers in {:.2?}",
        inputs.len(),
        all_minimizers.len(),
        total_time
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_creation() {
        let header = IndexHeader::new(31, 21);

        assert_eq!(header.format_version, 2);
        assert_eq!(header.kmer_length(), 31);
        assert_eq!(header.window_size(), 21);
    }

    #[test]
    fn test_header_validation() {
        // Valid header
        let valid_header = IndexHeader {
            format_version: 2,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(valid_header.validate().is_ok());

        // Invalid format version
        let invalid_header = IndexHeader {
            format_version: 1, // Unsupported version
            kmer_length: 31,
            window_size: 21,
        };
        assert!(invalid_header.validate().is_err());
    }
}
