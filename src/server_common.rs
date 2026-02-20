//! Common structures and types used in the client and server
use crate::IndexHeader;
use anyhow::Result;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};

/// Request structure for filtering sequences
#[derive(Serialize, Deserialize)]
pub struct FilterSequencesRequest {
    /// Sequences to filter
    pub sequences: Vec<Vec<u8>>,

    /// Absoulte filtering threshold
    pub abs_threshold: usize,

    /// Relative filtering threshold (proportion of minimizers that must match)
    pub rel_threshold: f64,

    /// Prefix length for minimizer computation
    pub prefix_length: usize,

    /// Whether running in deplete mode
    pub deplete: bool,

    /// Whether running in debug mode
    pub debug: bool,
}

/// Response structure for sequence filtering results
#[derive(Serialize, Deserialize)]
pub struct FilterSequencesResponse {
    /// Results for each sequence
    /// Each tuple contains:
    /// - A boolean indicating if the sequence should be output
    /// - If running in debug mode, a vector of strings of matched kmers. Empty vector else.
    /// - The length of the sequence
    pub results: Vec<(bool, Vec<String>, usize)>,
}

/// Request structure for filtering paired sequences
#[derive(Serialize, Deserialize)]
pub struct FilterPairedSequencesRequest {
    /// Sequences to filter
    pub sequences: Vec<(Vec<u8>, Vec<u8>)>,

    /// Absoulte filtering threshold
    pub abs_threshold: usize,

    /// Relative filtering threshold (proportion of minimizers that must match)
    pub rel_threshold: f64,

    /// Prefix length for minimizer computation
    pub prefix_length: usize,

    /// Whether running in deplete mode
    pub deplete: bool,

    /// Whether running in debug mode
    pub debug: bool,
}

/// Response structure for sequence filtering paired results
#[derive(Serialize, Deserialize)]
pub struct FilterPairedSequencesResponse {
    /// Results for each pair of sequences
    /// Each tuple contains:
    /// - A boolean indicating if the pair should be output
    /// - If running in debug mode, a vector of strings of matched kmers. Empty vector else.
    /// - The length of the first sequence
    /// - The length of the second sequence
    pub results: Vec<(bool, Vec<String>, usize, usize)>,
}

/// Get the header of the index loaded into a remote server
/// Required in order to ensure that the locally computed minimizers match
/// the kmer length and window size
pub fn get_server_index_header(server_address: &str) -> Result<IndexHeader> {
    // Create a client to send the minimizers to the server
    let client = Client::new();

    // Send the minimizers as a POST request
    let response = client
        .get(server_address.to_owned() + "/index_header")
        .send()?;

    // Check if the response indicates a match
    if response.status().is_success() {
        Ok(response.json::<IndexHeader>()?)
    } else {
        Err(anyhow::anyhow!(
            "Server returned an error: {}",
            response.status()
        ))
    }
}
